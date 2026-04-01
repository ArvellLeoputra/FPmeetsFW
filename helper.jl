# Helper function to check constraint feasibility
function check_feasibility(
    scip::Ptr{SCIP.SCIP_},
    solution::Vector{Float64},
    col_to_idx::Dict{Ptr{SCIP.SCIP_COL}, Int},
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)

    ptr_rows = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, ncols)

    # Check bounds
    for j in 1:ncols
        var = SCIP.SCIPcolGetVar(lp_cols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if solution[j] < lb - tolerance || solution[j] > ub + tolerance
            return false
        end
    end

    # Constraint check using rows
    for i in 1:nrows
        row = lp_rows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonz_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonz_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        activity = 0.0
        
        for k in 1:nnonz
            col = nonz_cols[k]
            idx = col_to_idx[col]
            activity += nonz_vals[k] * solution[idx]
        end
    
        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

        if lhs > -SCIP.SCIPinfinity(scip) && activity < lhs - tolerance
            return false
        end

        if rhs < SCIP.SCIPinfinity(scip) && activity > rhs + tolerance
            return false
        end
    end

    return true
end

# Helper function to check integrality
function check_integrality(
    solution::Vector{Float64}, 
    binary_indices::Vector{Int}, 
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    for i in binary_indices
        if abs(solution[i] - round(solution[i])) > tolerance
            return false
        end
    end
    return true
end