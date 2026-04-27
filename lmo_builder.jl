# LMO Builder from SCIP LP
# TODO: Build LMO independently of SCIP LP
# TODO: Multi heuristic calls
function build_lmo_from_scip_lp(scip::Ptr{SCIP.SCIP_}, nvars, nrows)
    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

    ptr_rows = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    # Build LMO model
    opt_model = SCIP.Optimizer()
    MOI.set(opt_model, MOI.RawOptimizerAttribute("presolving/maxrounds"), 0)
    MOI.set(opt_model, MOI.RawOptimizerAttribute("display/verblevel"), 0)
    x = MOI.add_variables(opt_model, nvars)

    # Add variable bounds
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)

        # TODO: Multi heuristic calls
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        # if lb <= -SCIP.SCIPinfinity(scip) || ub >= SCIP.SCIPinfinity(scip)                                                                                                                                                                
        #     println("  var $j: lb=$lb  ub=$ub  (UNBOUNDED)")                                                                                                                                                                              
        # end
        
        if lb > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(opt_model, x[j], MOI.GreaterThan(lb))
        end
        if ub < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(opt_model, x[j], MOI.LessThan(ub))
        end
    end

    # Add constraints
    # TODO: Multi heuristic calls
    col_to_idx = Dict(lp_cols[k] => k for k in 1:nvars)

    for i in 1:nrows
        row = lp_rows[i]
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzero_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzero_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        terms = [MOI.ScalarAffineTerm(nonzero_vals[k], x[col_to_idx[nonzero_cols[k]]]) for k in 1:nnonz]
        aff = MOI.ScalarAffineFunction(terms, 0.0)

        # SCIP stores LP rows as lhs <= ax + const <= rhs
        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        if lhs > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(opt_model, aff, MOI.GreaterThan(lhs))
        end
        if rhs < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(opt_model, aff, MOI.LessThan(rhs))
        end
    end

    # use_modify = false to set a new objective each iteration without modifying the model structure
    return FrankWolfe.MathOptLMO(opt_model, false)  # might be slower, but safer
end
