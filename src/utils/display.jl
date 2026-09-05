mutable struct PumpDisplayColumn
    name::String
    width::Int
    decimals::Int
end

mutable struct PumpDisplay
    column::Vector{PumpDisplayColumn}
end

function addColumn!(display::PumpDisplay, name::String, width::Int, decimals::Int = 0)
    push!(display.column, PumpDisplayColumn(name, width, decimals))
end

function printHeader!(display::PumpDisplay)
    for col in display.column
        print(rpad(col.name, col.width))
    end
    println()
end

function formatValue(val::Float64, decimals::Int)
    if isnan(val)
        return "NaN"
    end

    fixedStr = Printf.format(Printf.Format("%.$(decimals)f"), val)
    dotIdx = findfirst('.', fixedStr)
    intDigits = count(isdigit, fixedStr[1:dotIdx-1])
    if intDigits <= DEF_MAX_INT_DIGITS
        return fixedStr
    end

    # Integer part exceeds the fixed digit cap: fall back to scientific notation,
    # which has a fixed length regardless of magnitude.
    return Printf.format(Printf.Format("%.$(decimals)e"), val)
end

function formatFlips(flips::Int)
    return flips > 0 ? string(flips) : ""
end

function printRow!(display::PumpDisplay, values...)
    # Reprint the header every 10 iterations (values[1] is the iteration number).
    iter = values[1]
    if iter isa Integer && iter > 1 && (iter - 1) % 10 == 0
        printHeader!(display)
    end

    for (col, val) in zip(display.column, values)
        if val isa Float64
            val = formatValue(val, col.decimals)
        end
        print(rpad(string(val), col.width))
    end
    println()
end

function printInitialSolveInfo(scip, initObj, intIdx)
    # deg measures how degenerate the face is (fraction of non-basic vars with zero reduced cost);
    # high degeneracy => the projection LMO can jump vertices between pump iterations.
    # varconsratio estimates the face's dimension
    deg = Ref{Cdouble}(0.0)
    varconsratio = Ref{Cdouble}(0.0)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPDualDegeneracy(scip, deg, varconsratio)

    printstyled("[initialSolve]\n", color=:cyan)
    @printf("presolve = %d->%d vars, %d->%d conss, %d rounds, %.2fs\n",
        SCIP.SCIPgetNOrigVars(scip), SCIP.SCIPgetNVars(scip),
        SCIP.SCIPgetNOrigConss(scip), SCIP.SCIPgetNConss(scip),
        SCIP.SCIPgetNPresolRounds(scip), SCIP.SCIPgetPresolvingTime(scip))
    @printf("rootLP = %dx%d, origObj %.6f, %d iters, %.2fs\n",
        SCIP.SCIPgetNLPCols(scip), SCIP.SCIPgetNLPRows(scip),
        initObj, SCIP.SCIPgetNLPIterations(scip), SCIP.SCIPgetFirstLPTime(scip))
    @printf("lpSolves = %d (primal %d, dual %d, resolve %d)\n",
        SCIP.SCIPgetNLPs(scip), SCIP.SCIPgetNPrimalLPs(scip),
        SCIP.SCIPgetNDualLPs(scip), SCIP.SCIPgetNResolveLPs(scip))
    @printf("fracVars = %d / %d\n", SCIP.SCIPgetNLPBranchCands(scip), length(intIdx))
    @printf("separation = %d rounds, %d cuts\n",
        SCIP.SCIPgetNSepaRounds(scip), SCIP.SCIPgetNCutsApplied(scip))
    @printf("degeneracy = %.4f, var/cons %.4f\n", deg[], varconsratio[])
end

function printInitialBasisInfo(cstat, rstat)
    @printf("Initial basis: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d\n",
        count(==(0), cstat), count(==(1), cstat), count(==(2), cstat),
        count(==(0), rstat), count(==(1), rstat), count(==(2), rstat))
end

function setupPumpDisplay()
    pumpDisplay = PumpDisplay(PumpDisplayColumn[])
    addColumn!(pumpDisplay, "iter", 7)
    addColumn!(pumpDisplay, "stage", 7)
    addColumn!(pumpDisplay, "origObj", 15, 2)
    addColumn!(pumpDisplay, "projObj", 15, 4)
    addColumn!(pumpDisplay, "step", 15, 4)
    addColumn!(pumpDisplay, "nFrac", 8)
    addColumn!(pumpDisplay, "fwIters", 10)
    addColumn!(pumpDisplay, "time", 10, 2)
    addColumn!(pumpDisplay, "#flips", 8)
    addColumn!(pumpDisplay, "P", 3)
    addColumn!(pumpDisplay, "R", 3)
    addColumn!(pumpDisplay, "status", 14)

    printstyled("[pump]\n", color=:cyan)
    printHeader!(pumpDisplay)

    return pumpDisplay
end

# Verbose-mode counterpart to printRow! for a single FW-projection iteration
function printVerboseIteration(
    stage::Int,
    origObj::Float64,
    projObj::Float64,
    step::Float64,
    nFrac::Int,
    fwIters::Int,
    iterTime::Float64,
    flips::String,
    perturbed::Bool,
    restarted::Bool,
    outcome::String
)
    @printf("FW: stage=%d origObj=%.4f projObj=%.4f step=%.4f nFrac=%d fwIters=%d iterTime=%.4fs #flips=%s P=%s R=%s -> %s\n",
        stage, origObj, projObj, step, nFrac, fwIters, iterTime, flips,
        perturbed ? "*" : " ", restarted ? "*" : " ", outcome)
end

# Logs one FW-projection iteration's outcome
function logIteration(
    config, pumpDisplay, stats, stage, origObj, projObj, step, nFrac, fwIters,
    iterTime, heurStartTime, flips, perturbed, restarted, rowLabel, verboseLabel
)
    if config.verbose == 1
        printRow!(pumpDisplay, stats.pumpIterations, stage, origObj, projObj, step, nFrac, fwIters,
            timeElapsed(heurStartTime), formatFlips(flips), perturbed ? "*" : "", restarted ? "*" : "", rowLabel)
    
    elseif config.verbose >= 2
        printVerboseIteration(stage, origObj, projObj, step, nFrac, fwIters, iterTime,
            formatFlips(flips), perturbed, restarted, verboseLabel)
    end
end

# Logs a direct-acceptance event (RandFeasCheck / FeasRound / DiveSolve)
function logDirectAccept(config, pumpDisplay, stats, stage, origObj, step, heurStartTime, flips, perturbed, restarted, rowLabel, verboseLabel)
    elapsed = timeElapsed(heurStartTime)
    if config.verbose == 1
        printRow!(pumpDisplay, stats.pumpIterations, stage, origObj, 0.0, step, 0, 0, elapsed,
            formatFlips(flips), perturbed ? "*" : "", restarted ? "*" : "", rowLabel)

    elseif config.verbose >= 2
        @printf("%s: stage=%d origObj=%.4f step=%.4f elapsed=%.4fs #flips=%s P=%s R=%s\n",
            verboseLabel, stage, origObj, step, elapsed, formatFlips(flips), perturbed ? "*" : " ", restarted ? "*" : " ")
    end
end