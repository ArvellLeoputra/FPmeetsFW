function getLPInfo(scip::Ptr{SCIP.SCIP_})
    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    colsPtr = SCIP.SCIPgetLPCols(scip)
    rowsPtr = SCIP.SCIPgetLPRows(scip)

    lpCols = copy(unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, colsPtr, ncols))
    lpRows = copy(unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, rowsPtr, nrows))
    colDict = Dict(lpCols[k] => k for k in 1:ncols)

    binIdx = Int[]
    gIntIdx = Int[]
    objCoeffs = zeros(SCIP.SCIP_Real, ncols)
    initSol = zeros(SCIP.SCIP_Real, ncols)

    for j in 1:ncols
        var = SCIP.SCIPcolGetVar(lpCols[j])
        objCoeffs[j] = SCIP.SCIPvarGetObj(var)

        # Remove fixed variables from the index lists
        isFixed = SCIP.SCIPvarGetLbLocal(var) == SCIP.SCIPvarGetUbLocal(var)
        if !isFixed && SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_BINARY
            push!(binIdx, j)
        elseif !isFixed && SCIP.SCIPvarGetType(var) == SCIP.SCIP_VARTYPE_INTEGER
            push!(gIntIdx, j)
        end
        
        initSol[j] = SCIP.SCIPcolGetPrimsol(lpCols[j])
    end

    binSet = Set(binIdx)
    intIdx = [binIdx; gIntIdx]
    objScale = sqrt(sum(abs2, objCoeffs))

    return LPInfo(; lpCols, lpRows, objCoeffs, objScale, colDict, binIdx,
                  binSet, gIntIdx, intIdx, ncols, nrows, initSol)
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

    inf = SCIP.SCIPinfinity(scip)

    # Constraint check using rows
    for i in 1:length(lpRows)
        row = lpRows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        activity = 0.0
        for k in 1:nnonz
            activity += nonzVals[k] * sol[colDict[nonzCols[k]]]
        end

        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant
        
        belowLhs = lhs > -inf && SCIP.SCIPisFeasLT(scip, activity, lhs) == SCIP.TRUE
        aboveRhs = rhs < inf && SCIP.SCIPisFeasGT(scip, activity, rhs) == SCIP.TRUE

        if belowLhs || aboveRhs
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