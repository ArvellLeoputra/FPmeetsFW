function timeElapsed(start::Float64)
    return time() - start
end

function buildStats(scip::SCIP.SCIPData, heur::FPFWHeuristic, totalTime::Float64)
    stats = heur.data.stats
    if heur.data.called > 0
        stats.totalTime = totalTime
        stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)
        return stats
    end

    stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))
    stats.dualBound = Float64(SCIP.SCIPgetDualbound(scip))
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    stats.totalTime = totalTime
    stats.solutionFound = SCIP.SCIPgetNSols(scip) > 0

    status = SCIP.SCIPgetStatus(scip)
    if status == SCIP.SCIP_STATUS_OPTIMAL
        stats.exitReason = :scip_optimal
    elseif status == SCIP.SCIP_STATUS_INFEASIBLE
        stats.exitReason = :scip_infeasible
    elseif status == SCIP.SCIP_STATUS_UNBOUNDED
        stats.exitReason = :scip_unbounded
    elseif status == SCIP.SCIP_STATUS_TIMELIMIT
        stats.exitReason = :scip_time_limit
    elseif status == SCIP.SCIP_STATUS_NODELIMIT
        stats.exitReason = :scip_node_limit
    else
        stats.exitReason = :scip_unknown
    end

    stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)

    return stats
end

function printRunInfo(scip::SCIP.SCIPData)
    name = unsafe_string(SCIP.SCIPgetProbName(scip))
    varCount = SCIP.SCIPgetNOrigVars(scip)
    binCount = SCIP.SCIPgetNBinVars(scip)
    intCount = SCIP.SCIPgetNIntVars(scip)
    contCount = SCIP.SCIPgetNContVars(scip)

    printstyled("[run info]\n", color=:cyan)
    println("instance = $name")
    println("totalVars = $varCount")
    println("binaryVars = $binCount")
    println("integerVars = $intCount")
    println("continuousVars = $contCount")
end

function printConfigs(config::FPFWConfig)
    printstyled("[FPFW configs]\n", color=:cyan)
    println("norm = $(config.norm)")
    println("fwVariant = $(config.fwVariant)")
    println("lineSearch = $(config.lineSearch)")
    println("timeLimit = $(config.timeLimit)")
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
    println("fwWarmStart = $(config.fwWarmStart ? "enabled" : "disabled")")
    println("lmoWarmStart = $(config.lmoWarmStart ? "enabled" : "disabled")")
    println("useSubMIP = $(config.useSubMIP ? "enabled" : "disabled")")
    println("presolve = $(config.presolve ? "enabled" : "disabled")")
    println("seed = $(config.seed)")
end

function printResults(stats::FPFWStats)
    exitMsg = if stats.exitReason == :time_limit
        "time limit reached"
    elseif stats.exitReason == :infeasible_fw
        "FW returned a point outside the feasible polytope (numerical error)"
    elseif stats.exitReason == :solution_found
        "integer feasible solution accepted by SCIP at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :rr_solution_found
        "integer feasible solution found by randomized rounding at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :scip_optimal
        "problem solved to optimality by SCIP before heuristic was called"
    elseif stats.exitReason == :scip_infeasible
        "LP is infeasible, thus MIP is infeasible"
    elseif stats.exitReason == :scip_time_limit
        "time limit reached before heuristic was called"
    elseif stats.exitReason == :scip_node_limit
        "SCIP node limit reached before heuristic was called"
    elseif stats.exitReason == :scip_unbounded
        "LP is unbounded, thus MIP is unbounded"
    else
        "unknown exit"
    end

    printstyled("[result]\n", color=:cyan)
    println("primalBound = $(round(stats.primalBound, digits=4))")
    println("dualBound = $(round(stats.dualBound, digits=4))")
    println("gap = $(@sprintf("%.2f %%", stats.gap * 100))")
    println("primalIntegral = $(round(stats.primalIntegral, digits=4))")
    println("solFound = $(stats.solutionFound)")
    println("totalTime = $(round(stats.totalTime, digits=2))s")

    if !startswith(string(stats.exitReason), "scip_")
        println("totalHeurTime = $(round(stats.heurTime, digits=2))s")
        println("fwTime = $(round(stats.fwTime, digits=2))s")
        println("randRoundTime = $(round(stats.rrTime, digits=2))s")
        println("pumpIterations = $(stats.pumpIterations)")
        println("fwIterations = $(stats.fwIterations)")
        println("perturbCount = $(stats.perturbCount)")
        println("restartCount = $(stats.restartCount)")
    end

    println("exitReason = $exitMsg")
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

function printRow!(display::PumpDisplay, values...)
    for (col, val) in zip(display.column, values)
        if val isa Float64
            val = round(val, digits=col.decimals)
        end
        print(rpad(string(val), col.width))
    end
    println()
end

function areSolutionsEqual(scip::Ptr{SCIP.SCIP_}, intIdx::Vector{Int}, x1::Vector{Float64}, x2::Vector{Float64})
    for i in intIdx
        if SCIP.SCIPisEQ(scip, x1[i], x2[i]) == SCIP.FALSE
            return false
        end
    end
    return true
end

function countFracVars(scip::Ptr{SCIP.SCIP_}, intIdx::Vector{Int},x::Vector{Float64})
    cnt = 0
    for i in intIdx
        if SCIP.SCIPisEQ(scip, x[i], round(x[i])) == SCIP.FALSE
            cnt += 1
        end
    end
    return cnt
end

function getLPData(scip::Ptr{SCIP.SCIP_}, ncols::Int32, nrows::Int32)
    colsPtr = SCIP.SCIPgetLPCols(scip)
    rowsPtr = SCIP.SCIPgetLPRows(scip)
    
    lpCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, colsPtr, ncols)
    lpRows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, rowsPtr, nrows)
    colDict = Dict(lpCols[k] => k for k in 1:ncols)
    
    binIdx = Int[]
    intIdx = Int[]
    initSol = zeros(SCIP.SCIP_Real, ncols)

    for j in 1:ncols
        var = SCIP.SCIPcolGetVar(lpCols[j])
        if SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_BINARY
            push!(binIdx, j)
        elseif SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_INTEGER
            push!(intIdx, j)
        end
        initSol[j] = SCIP.SCIPcolGetPrimsol(lpCols[j])
    end

    return lpCols, lpRows, colDict, binIdx, intIdx, initSol
end

function computePrimalIntegral(events::Vector{Tuple{Float64, Float64}}, totalTime::Float64)
    area = 0.0
    prevTime = 0.0
    prevGap  = 1.0

    for (t, gap) in events
        area += prevGap * (t - prevTime)
        prevTime = t
        prevGap = gap
    end

    area += prevGap * (totalTime - prevTime)
    return area
end

# Rounding threshold generator
function getRoundingThreshold(random::Bool)
    return random ? rand() : 0.5
end

function roundSolution!(xRound::Vector{Float64}, x::Vector{Float64}, intIdx::Vector{Int}, randRound::Bool)
    threshold = getRoundingThreshold(randRound)
    for i in intIdx
        xRound[i] = floor(x[i] + threshold)
    end
end

# Hash function for cycle detection (only hashes integer variable values)
function hashRounded(x::Vector{Float64}, intIdx::Vector{Int})
    hash(tuple((x[i] for i in intIdx)...))
end

function perturb(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
    binIdx::Vector{Int},
    intIdx::Vector{Int},
    avgFlips::Int
)
    fracVars = [(abs(xRound[i] - x[i]), i) for i in intIdx if SCIP.SCIPisEQ(scip, xRound[i], x[i]) == SCIP.FALSE]
    sort!(fracVars, rev=true)

    nFlips = round(Int, avgFlips * (rand() + 0.5))
    nFracFlips = clamp(nFlips, 1, length(fracVars))

    if DEBUG_VERBOSE
        println("Perturbing $nFracFlips integer variables (out of $(length(fracVars)) fractional variables)")
    end

    binSet = Set(binIdx)  # faster lookup
    for k in 1:nFracFlips
        _, i = fracVars[k]
        if i in binSet  # binary: round to opposite
            xRound[i] = 1.0 - xRound[i]
        else  # general integer: reverse rounding direction
            xRound[i] += SCIP.SCIPisLT(scip, xRound[i], x[i]) == SCIP.TRUE ? 1.0 : -1.0
        end
    end
end

function restart(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
    prevRound::Vector{Float64},
    binIdx::Vector{Int},
    gIntIdx::Vector{Int},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    avgFlips::Int
)
    changed = 0

    # Binary variables
    for i in binIdx
        r = rand() - 0.47  # [-0.47, 0.53)

        if r > 0 && SCIP.SCIPisEQ(scip, xRound[i], prevRound[i]) == SCIP.TRUE  # stuck variable
            sigma = abs(xRound[i] - x[i])
            if sigma + r > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    if changed == 0
        for i in binIdx
            if rand() > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    # General integer variables
    if !isempty(gIntIdx)
        for _ in 1:avgFlips
            k = rand(1:length(gIntIdx))
            i = gIntIdx[k]

            var = SCIP.SCIPcolGetVar(lpCols[i])
            lb = SCIP.SCIPvarGetLbLocal(var)
            ub = SCIP.SCIPvarGetUbLocal(var)
            r = rand()

            newVal = if (ub - lb) < DEF_BIGBIGM
                floor(lb + (1 + ub - lb) * r)
            elseif (xRound[i] - lb) < DEF_BIGM
                lb + (2 * DEF_BIGM - 1) * r
            elseif (ub - xRound[i]) < DEF_BIGM
                ub - (2 * DEF_BIGM - 1) * r
            else
                xRound[i] + (2 * DEF_BIGM - 1) * r - DEF_BIGM
            end

            xRound[i] = clamp(floor(newVal), lb, ub)
            changed += 1
        end
    end
end

# Fix integer variables to xRound and solve the LP to adjust continuous variables
# Fallback when xRound alone is not accepted as a feasible MIP solution
function subMIPsolve(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    intIdx::Vector{Int},
    xRound::Vector{Float64},
    ncols::Int32
)::Tuple{Bool, Vector{Float64}}

    SCIP.SCIPstartDive(scip)
    try
        # Fix integer variables to their rounded values
        for i in intIdx
            var = SCIP.SCIPcolGetVar(lpCols[i])
            SCIP.SCIPchgVarLbDive(scip, var, xRound[i])
            SCIP.SCIPchgVarUbDive(scip, var, xRound[i])
        end

        lperror = Ref{SCIP.SCIP_Bool}(SCIP.FALSE)
        cutoff  = Ref{SCIP.SCIP_Bool}(SCIP.FALSE)

        SCIP.SCIPsolveDiveLP(scip, -1, lperror, cutoff)

        if lperror[] == SCIP.TRUE || cutoff[] == SCIP.TRUE
            return false, Float64[]
        end

        if SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL
            sol = [SCIP.SCIPcolGetPrimsol(lpCols[j]) for j in 1:ncols]
            return true, sol
        else
            return false, Float64[]
        end
    finally
        SCIP.SCIPendDive(scip)
    end
end

# Helper function to check LP feasibility
function isSolutionLPFeasible(
    scip::Ptr{SCIP.SCIP_},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    sol::Vector{Float64},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
)::Bool
    # Check bounds
    for j in 1:length(lpCols)
        var = SCIP.SCIPcolGetVar(lpCols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if SCIP.SCIPisFeasLT(scip, sol[j], lb) == SCIP.TRUE || SCIP.SCIPisFeasGT(scip, sol[j], ub) == SCIP.TRUE
            return false
        end
    end

    # Constraint check using rows
    for i in 1:length(lpRows)
        row = lpRows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        activity = 0.0
        
        for k in 1:nnonz
            col = nonzCols[k]
            idx = colDict[col]
            activity += nonzVals[k] * sol[idx]
        end
    
        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        if lhs > -SCIP.SCIPinfinity(scip) && SCIP.SCIPisFeasLT(scip, activity, lhs) == SCIP.TRUE
            return false
        end

        if rhs < SCIP.SCIPinfinity(scip) && SCIP.SCIPisFeasGT(scip, activity, rhs) == SCIP.TRUE
            return false
        end
    end

    return true
end

# Helper function to check integrality
function isSolutionIntegral(scip::Ptr{SCIP.SCIP_}, sol::Vector{Float64}, intIdx::Vector{Int})
    for i in intIdx
        if SCIP.SCIPisEQ(scip, sol[i], round(sol[i])) == SCIP.FALSE
            return false
        end
    end
    return true
end

function submitSolution(
    scip::Ptr{SCIP.SCIP_},
    heurPtr::Ptr{SCIP.SCIP_HEUR},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    sol::Vector{Float64},
    ncols::Int32
)
    solPtr = Ref{Ptr{SCIP.SCIP_SOL}}()
    SCIP.SCIPcreateSol(scip, solPtr, heurPtr)
    scipSol = solPtr[]

    # Set solution values
    for j in 1:ncols
        col = lpCols[j]
        var = SCIP.SCIPcolGetVar(col)
        SCIP.SCIPsetSolVal(scip, scipSol, var, sol[j])
    end

    # Try to add solution
    stored = Ref{SCIP.SCIP_Bool}()
    SCIP.SCIPtrySol(scip, scipSol, SCIP.FALSE, SCIP.FALSE, SCIP.TRUE, SCIP.TRUE, SCIP.TRUE, stored)

    if stored[] == SCIP.TRUE
        return true
    else
        SCIP.SCIPfreeSol(scip, solPtr)
        return false
    end
end