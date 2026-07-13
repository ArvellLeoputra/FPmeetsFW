# LMO Builder from SCIP LP
# TODO: Build LMO independently of SCIP LP
# TODO: Multi heuristic calls
function SCIPbuildLMO(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    ncols::Int32,
    nrows::Int32
)
    lmoStart = time()

    # Build LMO model
    optModel = SCIP.Optimizer()
    MOI.set(optModel, MOI.RawOptimizerAttribute("presolving/maxrounds"), 0)
    MOI.set(optModel, MOI.RawOptimizerAttribute("display/verblevel"), 0)
    x = MOI.add_variables(optModel, ncols)

    # Add variable bounds
    for j in 1:ncols
        col = lpCols[j]
        var = SCIP.SCIPcolGetVar(col)

        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if lb > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(optModel, x[j], MOI.GreaterThan(lb))
        end
        if ub < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(optModel, x[j], MOI.LessThan(ub))
        end
    end

    # Add constraints
    for i in 1:nrows
        row = lpRows[i]
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        terms = [MOI.ScalarAffineTerm(nonzVals[k], x[colDict[nonzCols[k]]]) for k in 1:nnonz]
        aff = MOI.ScalarAffineFunction(terms, 0.0)

        # SCIP stores LP rows as lhs <= ax + const <= rhs
        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        if lhs > -SCIP.SCIPinfinity(scip)
            MOI.add_constraint(optModel, aff, MOI.GreaterThan(lhs))
        end
        if rhs < SCIP.SCIPinfinity(scip)
            MOI.add_constraint(optModel, aff, MOI.LessThan(rhs))
        end
    end

    lmoTime = timeElapsed(lmoStart)
    println("LMO build time : $lmoTime")

    # use_modify = false to set a new objective each iteration without modifying the model structure
    return FrankWolfe.MathOptLMO(optModel, false)  # might be slower, but safer
end