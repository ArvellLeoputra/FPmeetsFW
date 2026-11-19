# MOI-based LMO (lmoWarmStart = false): a fresh SCIP model built from SCIP's LP, no warm start
# TODO: Build LMO independently of SCIP LP
# TODO: Multi heuristic calls
function SCIPbuildLMO(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    gIntIdx::Vector{Int},
    norm::Symbol,
    ncols::Int32,
    nrows::Int32
)
    optModel = SCIP.Optimizer()
    # Turn off presolve, output, and symmetry detection
    MOI.set(optModel, MOI.RawOptimizerAttribute("presolving/maxrounds"), 0)
    MOI.set(optModel, MOI.RawOptimizerAttribute("display/verblevel"), 0)
    MOI.set(optModel, MOI.RawOptimizerAttribute("misc/usesymmetry"), 0)
    inf = SCIP.SCIPinfinity(scip)
    x = MOI.add_variables(optModel, ncols)

    # Copy the variable bounds
    for j in 1:ncols
        var = SCIP.SCIPcolGetVar(lpCols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if lb > -inf
            MOI.add_constraint(optModel, x[j], MOI.GreaterThan(lb))
        end
        if ub < inf
            MOI.add_constraint(optModel, x[j], MOI.LessThan(ub))
        end
    end

    # Copy the rows, moving each row's constant into its sides
    for i in 1:nrows
        row = lpRows[i]
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)
        terms = [MOI.ScalarAffineTerm(nonzVals[k], x[colDict[nonzCols[k]]]) for k in 1:nnonz]
        aff = MOI.ScalarAffineFunction(terms, 0.0)

        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        if lhs > -inf
            MOI.add_constraint(optModel, aff, MOI.GreaterThan(lhs))
        end
        if rhs < inf
            MOI.add_constraint(optModel, aff, MOI.LessThan(rhs))
        end
    end

    auxConstraints = Tuple{MOI.ConstraintIndex, MOI.ConstraintIndex}[]
    if norm == :manhattan
        aux = MOI.add_variables(optModel, length(gIntIdx))
        for (k, i) in enumerate(gIntIdx)
            MOI.add_constraint(optModel, aux[k], MOI.GreaterThan(0.0))
            # Row 1: aux_k - x_i >= -xRound_i, i.e. aux_k >= x_i - xRound_i (xRound set in MOIupdateRounding!)
            c1 = MOI.add_constraint(
                optModel,
                MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, aux[k]), MOI.ScalarAffineTerm(-1.0, x[i])], 0.0),
                MOI.GreaterThan(0.0)
            )
            # Row 2: aux_k + x_i >= xRound_i, i.e. aux_k >= xRound_i - x_i
            c2 = MOI.add_constraint(
                optModel,
                MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, aux[k]), MOI.ScalarAffineTerm(1.0, x[i])], 0.0),
                MOI.GreaterThan(0.0)
            )
            push!(auxConstraints, (c1, c2))
        end
    end

    # use_modify = false: replace the whole objective on each solve
    return MOILMO(FrankWolfe.MathOptLMO(optModel, false), Ref(Inf)), auxConstraints
end

# Update the aux rows' sides for a new xRound
function MOIupdateRounding!(
    lmo::MOILMO,
    auxConstraints::Vector{Tuple{MOI.ConstraintIndex, MOI.ConstraintIndex}},
    gIntIdx::Vector{Int},
    xRound::Vector{Float64}
)
    for (k, i) in enumerate(gIntIdx)
        c1, c2 = auxConstraints[k]
        MOI.set(lmo.inner.o, MOI.ConstraintSet(), c1, MOI.GreaterThan(-xRound[i]))
        MOI.set(lmo.inner.o, MOI.ConstraintSet(), c2, MOI.GreaterThan(xRound[i]))
    end
end

# One LP solve, bounded by the deadline. Throws the same exceptions as the LPI oracle.
function FrankWolfe.compute_extreme_point(lmo::MOILMO, direction::AbstractVector; kwargs...)
    timeLeft = lmo.deadline[] - time()
    timeLeft <= 0 && throw(LMODeadlineExceeded())
    MOI.set(lmo.inner.o, MOI.TimeLimitSec(), timeLeft)

    v = try
        FrankWolfe.compute_extreme_point(lmo.inner, direction; kwargs...)
    catch e
        e isa InterruptException && rethrow()
        nothing  # no solution to read, the status below says why
    end

    status = MOI.get(lmo.inner.o, MOI.TerminationStatus())
    status == MOI.TIME_LIMIT && throw(LMODeadlineExceeded())
    (v === nothing || status ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)) && throw(LMOSolverError())
    return v
end
