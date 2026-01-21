# LMO Builder from SCIP LP
# TODO: Build LMO independently of SCIP LP
function build_lmo_from_scip_lp(scip::Ptr{SCIP.SCIP_}, nvars, nrows)
    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

    ptr_rows = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    # Build MOI model from current LP
    moi_model = MOI.Utilities.Model{Float64}()
    x = MOI.add_variables(moi_model, nvars)

    # Add variable bounds
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)

        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if lb > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(moi_model, x[j], MOI.GreaterThan(lb))
        end
        if ub < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(moi_model, x[j], MOI.LessThan(ub))
        end
    end

    # Add constraints
    col_to_idx = Dict(lp_cols[k] => k for k in 1:nvars)

    for i in 1:nrows
        row = lp_rows[i]
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzero_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzero_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        terms = [MOI.ScalarAffineTerm(nonzero_vals[k], x[col_to_idx[nonzero_cols[k]]]) for k in 1:nnonz]
        aff = MOI.ScalarAffineFunction(terms, 0.0)

        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

        if lhs > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(moi_model, aff, MOI.GreaterThan(lhs))
        end
        if rhs < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(moi_model, aff, MOI.LessThan(rhs))
        end
    end

    # Create LMO
    opt_model = GLPK.Optimizer()
    MOI.copy_to(opt_model, moi_model)
    return FrankWolfe.MathOptLMO(opt_model)
end
