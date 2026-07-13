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

function origObjective(scip::Ptr{SCIP.SCIP_}, lpCols::Vector{Ptr{SCIP.SCIP_COL}}, sol::Vector{Float64}, ncols::Int32)
    transObj = sum(sol[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
    return SCIP.SCIPretransformObj(scip, transObj)
end

function isVarInteger(scip::Ptr{SCIP.SCIP_}, v::Float64)
    return SCIP.SCIPisEQ(scip, v, round(v)) == SCIP.TRUE
end

function countFracVars(scip::Ptr{SCIP.SCIP_}, intIdx::Vector{Int},x::Vector{Float64})
    cnt = 0
    for i in intIdx
        if !isVarInteger(scip, x[i])
            cnt += 1
        end
    end
    return cnt
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
        if !isVarInteger(scip, sol[i])
            return false
        end
    end
    return true
end

# Currently unused
function areSolutionsEqual(scip::Ptr{SCIP.SCIP_}, intIdx::Vector{Int}, x1::Vector{Float64}, x2::Vector{Float64})
    for i in intIdx
        if SCIP.SCIPisEQ(scip, x1[i], x2[i]) == SCIP.FALSE
            return false
        end
    end
    return true
end