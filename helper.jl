function printHeurSummary(stats::FPFWStats)
    exit_msg = if stats.exitReason == :time_limit
        "global time limit $(DEF_GLOBAL_TIME_LIMIT)s reached"
    elseif stats.exitReason == :restart_limit
        "FP cycled $(DEF_MAX_RESTARTS) times without progress"
    elseif stats.exitReason == :infeasible_fw
        "FW returned a point outside the feasible polytope (numerical error)"
    elseif stats.exitReason == :solutionFound
        "integer feasible solution accepted by SCIP at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :rr_solution_found
        "integer feasible solution found by randomized rounding at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :solution_rejected
        "integer feasible solution found but rejected by SCIP"
    elseif stats.exitReason == :scip_presolved
        "problem solved by SCIP during presolve"
    else
        "unknown exit"
    end

    heurTime = stats.heurTime == 0.0 ? "N/A" : "$(round(stats.heurTime, digits=2))s"  # 0.0 means heuristic never ran
    primalBound = stats.primalBound === nothing ? "N/A" : "$(round(stats.primalBound, digits=4))"
    gap = isinf(stats.gap) || stats.gap > 1e15 ? "Infinite" : @sprintf("%.2f %%", stats.gap * 100)

    printstyled("[result]\n", color=:cyan)
    println("primalBound = $primalBound")
    println("dualBound = $(stats.dualBound)")
    println("gap = $gap")
    println("totalTime = $(round(stats.totalTime, digits=2))s")
    println("totalHeurTime = $heurTime")
    println("fwTime = $(round(stats.fwTime, digits=2))s")
    println("randRoundTime = $(round(stats.rrTime, digits=2))s")
    println("pumpIterations = $(stats.pumpIterations)")
    println("fwIterations = $(stats.fwIterations)")
    println("restartCount = $(stats.restartCount)")
    println("solFound = $(stats.solutionFound)")
    println("exitReason = $exit_msg")
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

function areValsEqual(x1::Float64, x2::Float64, tol::Float64=DEF_INT_TOLERANCE)
    return abs(x1 - x2) <= tol
end

function isLowerThan(x1::Float64, x2::Float64, tol::Float64=DEF_INT_TOLERANCE)
    return x1 < x2 - tol
end

function areSolutionsEqual(intIdx::Vector{Int}, x1::Vector{Float64}, x2::Vector{Float64}, tol::Float64=DEF_INT_TOLERANCE)
    for i in intIdx
        if !areValsEqual(x1[i], x2[i], tol)
            return false
        end
    end
    return true
end

function countFracVars(intIdx::Vector{Int},x::Vector{Float64}, tol::Float64=DEF_INT_TOLERANCE)
    cnt = 0
    for i in intIdx
        if !areValsEqual(x[i], round(x[i]), tol)
            cnt += 1
        end
    end
    return cnt
end

function getLPData(scip::Ptr{SCIP.SCIP_}, lpData::LPData)
    lpData.nvars = SCIP.SCIPgetNLPCols(scip)
    lpData.nrows = SCIP.SCIPgetNLPRows(scip)

    colsPtr = SCIP.SCIPgetLPCols(scip)
    lpData.lpCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, colsPtr, lpData.nvars)
    lpData.colDict = Dict(lpData.lpCols[k] => k for k in 1:lpData.nvars)

    rowsPtr = SCIP.SCIPgetLPRows(scip)
    lpData.lpRows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, rowsPtr, lpData.nrows)

    lpData.initSol = zeros(SCIP.SCIP_Real, lpData.nvars)

    for j in 1:lpData.nvars
        var = SCIP.SCIPcolGetVar(lpData.lpCols[j])
        if SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_BINARY
            push!(lpData.binIdx, j)
        elseif SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_INTEGER
            push!(lpData.gIntIdx, j)
        end
        lpData.initSol[j] = SCIP.SCIPcolGetPrimsol(lpData.lpCols[j])
    end

    lpData.intIdx = [lpData.binIdx; lpData.gIntIdx]
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

function lpDiving!(scip::Ptr{SCIP.SCIP_}, lpData::LPData, xRound::Vector{Float64})::Tuple{Bool, Vector{Float64}}
    lpModel = Model(SCIP.Optimizer)
    set_silent(lpModel)

    vars = Vector{VariableRef}(undef, lpData.nvars)
    for j in 1:lpData.nvars
        var = SCIP.SCIPcolGetVar(lpData.lpCols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)
        vars[j] = @variable(lpModel, lower_bound=lb, upper_bound=ub)
    end

    for i in lpData.intIdx
        set_lower_bound(vars[i], xRound[i])
        set_upper_bound(vars[i], xRound[i])
    end

    for row in lpData.lpRows
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        activity = sum(nonzVals[k] * vars[lpData.colDict[nonzCols[k]]] for k in 1:nnonz)

        if lhs > -SCIP.SCIPinfinity(scip) && rhs < SCIP.SCIPinfinity(scip)
            @constraint(lpModel, lhs <= activity <= rhs)
        elseif lhs > -SCIP.SCIPinfinity(scip)
            @constraint(lpModel, activity >= lhs)
        elseif rhs < SCIP.SCIPinfinity(scip)
            @constraint(lpModel, activity <= rhs)
        end
    end

    @objective(lpModel, Min, sum(SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpData.lpCols[j])) * vars[j] for j in 1:lpData.nvars))

    optimize!(lpModel)

    if termination_status(lpModel) == MOI.OPTIMAL
        sol = [value(vars[j]) for j in 1:lpData.nvars]
        return true, sol
    else
        return false, Float64[]
    end
end

# Helper function to check LP feasibility
function isSolutionLPFeasible(scip::Ptr{SCIP.SCIP_}, lpData::LPData, sol::Vector{Float64}, tol::Float64=DEF_INT_TOLERANCE)
    # Check bounds
    for j in 1:length(lpData.lpCols)
        var = SCIP.SCIPcolGetVar(lpData.lpCols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if sol[j] < lb - tol || sol[j] > ub + tol
            return false
        end
    end

    # Constraint check using rows
    for i in 1:length(lpData.lpRows)
        row = lpData.lpRows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        activity = 0.0
        
        for k in 1:nnonz
            col = nonzCols[k]
            idx = lpData.colDict[col]
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
        if !areValsEqual(sol[i], round(sol[i]), tol)
            return false
        end
    end
    return true
end

function submitSolution(scip::Ptr{SCIP.SCIP_}, lpData::LPData, sol::Vector{Float64})
    return isSolutionIntegral(sol, lpData.intIdx) && isSolutionLPFeasible(scip, lpData, sol)
end