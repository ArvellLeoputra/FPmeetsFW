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

        # TODO: Multi heuristic calls
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
    # TODO: Multi heuristic calls
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
    println("LMO build time = $lmoTime")

    # use_modify = false to set a new objective each iteration without modifying the model structure
    return FrankWolfe.MathOptLMO(optModel, false)  # might be slower, but safer
end

function buildLPILMO(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    ncols::Int32,
    nrows::Int32
)
    lmoStart = time()

    # Create SCIPlpi instance
    lpiPtr = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    msghdlr = SCIP.SCIPgetMessagehdlr(scip)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpiPtr, msghdlr, "lmo", SCIP.SCIP_OBJSEN_MINIMIZE)
    lpi = lpiPtr[]

    obj = zeros(Cdouble, ncols)
    lb = zeros(Cdouble, ncols)
    ub = zeros(Cdouble, ncols)
    inf = SCIP.SCIPlpiInfinity(lpi)

    for j in 1:ncols
        col = lpCols[j]
        var = SCIP.SCIPcolGetVar(col)
        lb[j] = max(SCIP.SCIPvarGetLbLocal(var), -inf)
        ub[j] = min(SCIP.SCIPvarGetUbLocal(var), inf)
    end

    # Add columns
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(
        lpi,
        Cint(ncols),
        obj,
        lb,
        ub,
        C_NULL,
        Cint(0),
        C_NULL,
        C_NULL,
        C_NULL
    )

    lhsVec = Float64[]
    rhsVec = Float64[]
    begVec = Cint[]
    nonzIdx = Cint[]
    nonzVals = Float64[]
    offset = Cint(0)

    for i in 1:nrows
        row = lpRows[i]
        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant
        push!(lhsVec, lhs < -inf ? -inf : lhs)
        push!(rhsVec, rhs > inf ? inf : rhs)
        push!(begVec, offset)

        nnonz = SCIP.SCIProwGetNNonz(row)
        if nnonz > 0    
            rowCols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
            rowVals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)
            for k in 1:nnonz
                push!(nonzIdx, Cint(colDict[rowCols[k]] - 1))
                push!(nonzVals, rowVals[k])
            end
            offset += nnonz
        end
    end

    # Add rows
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(
        lpi,
        nrows,
        lhsVec,
        rhsVec,
        C_NULL,
        offset,
        begVec,
        nonzIdx,
        nonzVals
    )

    lmoTime = timeElapsed(lmoStart)
    println("LMO build time = $lmoTime")

    ncols_ref = Ref{Cint}(0)
    nrows_ref = Ref{Cint}(0)
    SCIP.SCIPlpiGetNCols(lpi, ncols_ref)
    SCIP.SCIPlpiGetNRows(lpi, nrows_ref)
    @assert ncols_ref[] == ncols "LPI ncols $(ncols_ref[]) != ncols $ncols"
    @assert nrows_ref[] == nrows "LPI nrows $(nrows_ref[]) != nrows $nrows"

    return BaseLMO(lpi, ncols, nrows)
end

# Extract LP basis from SCIP's internal LP solver
function LPIgetBase(scip::Ptr{SCIP.SCIP_}, ncols::Int32, nrows::Int32)
    # Get SCIP's internal LP handle
    lpiPtr = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPI(scip, lpiPtr)
    lpi = lpiPtr[]

    # Sanity check
    ncols_ref = Ref{Cint}(0)
    nrows_ref = Ref{Cint}(0)
    SCIP.SCIPlpiGetNCols(lpi, ncols_ref)
    SCIP.SCIPlpiGetNRows(lpi, nrows_ref)
    @assert ncols_ref[] == ncols "LPI ncols $(ncols_ref[]) != ncols $ncols"
    @assert nrows_ref[] == nrows "LPI nrows $(nrows_ref[]) != nrows $nrows"

    isOptimal = SCIP.SCIPlpiIsOptimal(lpi)
    @assert isOptimal == SCIP.TRUE "LPI is not optimal"

    # Extract current LP basis
    cstat = zeros(Cint, ncols)
    rstat = zeros(Cint, nrows)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetBase(lpi, cstat, rstat)

    return cstat, rstat
end

# Set LP basis
function LPIsetBase(lmo::BaseLMO, cstat::Vector{Cint}, rstat::Vector{Cint})
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetBase(lmo.lpi, cstat, rstat)
end

# Initialize LMO basis from SCIP's internal LP basis
function LPIinitBase(scip::Ptr{SCIP.SCIP_}, lmo::BaseLMO)
    cstat, rstat = LPIgetBase(scip, lmo.ncols, lmo.nrows)
    LPIsetBase(lmo, cstat, rstat)
end

function FrankWolfe.compute_extreme_point(lmo::BaseLMO, direction::AbstractVector; kwargs...)
    # Set objective
    SCIP.@SCIP_CALL SCIP.SCIPlpiChgObj(
        lmo.lpi,
        Cint(lmo.ncols),
        [Cint(i) for i in 0:(lmo.ncols - 1)],
        Vector{Cdouble}(direction)
    )

    # Check current basis statuses
    if DEBUG_VERBOSE
        nrows_ref = Ref{Cint}(0)
        SCIP.SCIPlpiGetNRows(lmo.lpi, nrows_ref)

        cstat = zeros(Cint, lmo.ncols)
        rstat = zeros(Cint, nrows_ref[])
        SCIP.SCIPlpiGetBase(lmo.lpi, cstat, rstat)

        @printf("LMO basis before solve: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d\n",
            count(==(0), cstat), count(==(1), cstat), count(==(2), cstat),
            count(==(0), rstat), count(==(1), rstat), count(==(2), rstat))
    end

    # Solve with dual simplex
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolveDual(lmo.lpi)

    if DEBUG_VERBOSE
        simplexIter = Ref{Cint}(0)
        SCIP.@SCIP_CALL SCIP.SCIPlpiGetIterations(lmo.lpi, simplexIter)

        nrows_ref = Ref{Cint}(0)
        SCIP.SCIPlpiGetNRows(lmo.lpi, nrows_ref)
        cstat_post = zeros(Cint, lmo.ncols)
        rstat_post = zeros(Cint, nrows_ref[])
        SCIP.SCIPlpiGetBase(lmo.lpi, cstat_post, rstat_post)

        @printf("LMO basis after solve: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d iters=%d\n",
            count(==(0), cstat_post), count(==(1), cstat_post), count(==(2), cstat_post),
            count(==(0), rstat_post), count(==(1), rstat_post), count(==(2), rstat_post),
            simplexIter[])
    end

    if SCIP.SCIPlpiIsPrimalFeasible(lmo.lpi) == SCIP.TRUE
        solVector = zeros(SCIP.SCIP_Real, lmo.ncols)
        obj = Ref{SCIP.SCIP_Real}(0)
        SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(
            lmo.lpi,
            obj,
            solVector,
            C_NULL,
            C_NULL,
            C_NULL
        )

        return solVector
    else
        error("BaseLMO: LP not primal feasible after dual simplex solve")
    end
end