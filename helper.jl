function buildStats(scip::SCIP.SCIPData, startTime::Float64, heur::FPFWHeuristic)
    if heur.called > 0
        return heur.stats
    end

    stats = FPFWStats()
    stats.primalBound = normalizeInf(Float64(SCIP.SCIPgetPrimalbound(scip)))
    stats.dualBound = normalizeInf(Float64(SCIP.SCIPgetDualbound(scip)))
    stats.gap = normalizeInf(Float64(SCIP.SCIPgetGap(scip)))

    stats.totalTime = time() - startTime
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

    return stats
end

function printRunInfo(scip::SCIP.SCIPData, fileName::String)
    varCount = SCIP.SCIPgetNOrigVars(scip)
    origVars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetOrigVars(scip), varCount)
    binCount = sum(SCIP.SCIPvarGetType(origVars[j]) == SCIP.SCIP_VARTYPE_BINARY for j in 1:varCount)
    intCount = sum(SCIP.SCIPvarGetType(origVars[j]) == SCIP.SCIP_VARTYPE_INTEGER for j in 1:varCount)
    contCount = varCount - binCount - intCount

    name = basename(fileName)
    while endswith(name, ".mps") || endswith(name, ".gz")
        name = splitext(name)[1]
    end

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
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
    println("warmStart = $(config.warmStart ? "enabled" : "disabled")")
    println("presolve = $(config.presolve ? "enabled" : "disabled")")
    println("seed = $(config.seed)")
end

function printResults(stats::FPFWStats)
    exitMsg = if stats.exitReason == :time_limit
        "global time limit $(DEF_GLOBAL_TIME_LIMIT)s reached"
    elseif stats.exitReason == :restart_limit
        "FP cycled $(DEF_MAX_RESTARTS) times without progress"
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
        "SCIP time limit $(DEF_SCIP_TIME_LIMIT)s reached before heuristic was called"
    elseif stats.exitReason == :scip_node_limit
        "SCIP node limit reached before heuristic was called"
    elseif stats.exitReason == :scip_unbounded
        "LP is unbounded, thus MIP is unbounded"
    else
        "unknown exit"
    end

    primalBound = isinf(stats.primalBound) ? "Inf" : "$(round(stats.primalBound, digits=4))"
    dualBound = isinf(stats.dualBound) ? "Inf" : "$(round(stats.dualBound, digits=4))"
    gap = isinf(stats.gap) ? "Inf" : @sprintf("%.2f %%", stats.gap * 100)

    printstyled("[result]\n", color=:cyan)
    println("primalBound = $primalBound")
    println("dualBound = $dualBound")
    println("gap = $gap")
    println("totalTime = $(round(stats.totalTime, digits=2))s")
    println("solFound = $(stats.solutionFound)")

    if !startswith(string(stats.exitReason), "scip_")
        println("totalHeurTime = $(round(stats.heurTime, digits=2))s")
        println("fwTime = $(round(stats.fwTime, digits=2))s")
        println("randRoundTime = $(round(stats.rrTime, digits=2))s")
        println("pumpIterations = $(stats.pumpIterations)")
        println("fwIterations = $(stats.fwIterations)")
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

function areValuesEqual(x1::Float64, x2::Float64, tol::Float64=DEF_INT_TOLERANCE)
    return abs(x1 - x2) <= tol
end

function isLowerThan(x1::Float64, x2::Float64, tol::Float64=DEF_INT_TOLERANCE)
    return x1 < x2 - tol
end

function areSolutionsEqual(intIdx::Vector{Int}, x1::Vector{Float64}, x2::Vector{Float64}, tol::Float64=DEF_INT_TOLERANCE)
    for i in intIdx
        if !areValuesEqual(x1[i], x2[i], tol)
            return false
        end
    end
    return true
end

function normalizeInf(x::Float64)
    if x > 1e15
        return Inf
    elseif x < -1e15
        return -Inf
    else
        return x
    end
end

function countFracVars(intIdx::Vector{Int},x::Vector{Float64}, tol::Float64=DEF_INT_TOLERANCE)
    cnt = 0
    for i in intIdx
        if !areValuesEqual(x[i], round(x[i]), tol)
            cnt += 1
        end
    end
    return cnt
end

function getLPData(scip::Ptr{SCIP.SCIP_}, nvars::Int32, nrows::Int32)
    colsPtr = SCIP.SCIPgetLPCols(scip)
    lpCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, colsPtr, nvars)
    colDict = Dict(lpCols[k] => k for k in 1:nvars)

    rowsPtr = SCIP.SCIPgetLPRows(scip)
    lpRows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, rowsPtr, nrows)

    binIdx = Int[]
    intIdx = Int[]
    initSol = zeros(SCIP.SCIP_Real, nvars)

    for j in 1:nvars
        var = SCIP.SCIPcolGetVar(lpCols[j])
        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binIdx, j)
        elseif SCIP.SCIPvarIsIntegral(var) == SCIP.TRUE
            push!(intIdx, j)
        end
        initSol[j] = SCIP.SCIPcolGetPrimsol(lpCols[j])
    end

    return lpCols, lpRows, colDict, binIdx, intIdx, initSol
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
function hashSolution(x::Vector{Float64}, intIdx::Vector{Int})
    hash(tuple((x[i] for i in intIdx)...))
end

function perturbSolution!(x::Vector{Float64}, xRound::Vector{Float64}, binIdx::Vector{Int}, gIntIdx::Vector{Int}, lpCols::Vector{Ptr{SCIP.SCIP_COL}})
    for i in binIdx
        if rand() < DEF_PERTURB_FRACTION
            xRound[i] = 1.0 - xRound[i]
        end
    end

    for i in gIntIdx
        if rand() < DEF_PERTURB_FRACTION
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
                x[i] + (2 * DEF_BIGM - 1) * r - DEF_BIGM
            end

            xRound[i] = clamp(floor(newVal), lb, ub)
        end
    end
end

function lpDiving!(scip::Ptr{SCIP.SCIP_}, lpCols::Vector{Ptr{SCIP.SCIP_COL}}, intIdx::Vector{Int}, xRound::Vector{Float64}, nvars::Int32)::Tuple{Bool, Vector{Float64}}
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
            sol = [SCIP.SCIPcolGetPrimsol(lpCols[j]) for j in 1:nvars]
            return true, sol
        else
            return false, Float64[]
        end
    finally
        SCIP.SCIPendDive(scip)
    end
end

# Helper function to check LP feasibility
function isSolutionLPFeasible(scip::Ptr{SCIP.SCIP_}, lpRows::Vector{Ptr{SCIP.SCIP_ROW}}, lpCols::Vector{Ptr{SCIP.SCIP_COL}}, sol::Vector{Float64}, colDict::Dict{Ptr{SCIP.SCIP_COL}, Int}, tol::Float64=DEF_INT_TOLERANCE)
    # Check bounds
    for j in 1:length(lpCols)
        var = SCIP.SCIPcolGetVar(lpCols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if sol[j] < lb - tol || sol[j] > ub + tol
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

        if lhs > -SCIP.SCIPinfinity(scip) && activity < lhs - tol
            return false
        end

        if rhs < SCIP.SCIPinfinity(scip) && activity > rhs + tol
            return false
        end
    end

    return true
end

# Helper function to check integrality
function isSolutionIntegral(sol::Vector{Float64}, intIdx::Vector{Int}, tol::Float64=DEF_INT_TOLERANCE)
    for i in intIdx
        if !areValuesEqual(sol[i], round(sol[i]), tol)
            return false
        end
    end
    return true
end

function submitSolution(scip::Ptr{SCIP.SCIP_}, heurPtr::Ptr{SCIP.SCIP_HEUR}, lpCols::Vector{Ptr{SCIP.SCIP_COL}}, sol::Vector{Float64}, nvars::Int32)
    solPtr = Ref{Ptr{SCIP.SCIP_SOL}}()
    SCIP.SCIPcreateSol(scip, solPtr, heurPtr)
    scipSol = solPtr[]

    # Set solution values
    for j in 1:nvars
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