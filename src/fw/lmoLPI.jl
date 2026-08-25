# The LMO's polytope (lmo.lpi) is built once per stage in buildLPILMO with a placeholder objective
function buildLPILMO(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    gIntIdx::Vector{Int},
    norm::Symbol,  # only :manhattan gets the auxiliary block
    stage::Int,
    ncols::Int32,
    nrows::Int32,
    verbose::Int
)
    # Create SCIPlpi instance
    lpiPtr = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    msghdlr = SCIP.SCIPgetMessagehdlr(scip)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpiPtr, msghdlr, "lmo", SCIP.SCIP_OBJSEN_MINIMIZE)
    lpi = lpiPtr[]

    obj = zeros(Cdouble, ncols)  # later set to the current gradient direction in compute_extreme_point
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

    # Add auxiliary block columns for Manhattan norm
    if norm == :manhattan
        nGInt = length(gIntIdx)
        auxObj = ones(Cdouble, nGInt)
        auxLb = zeros(Cdouble, nGInt)
        auxUb = fill(SCIP.SCIPlpiInfinity(lpi), nGInt)
        
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(
            lpi,
            Cint(nGInt),
            auxObj,
            auxLb,
            auxUb,
            C_NULL,
            Cint(0),
            C_NULL,
            C_NULL,
            C_NULL
        )

        auxLhs = Float64[]
        auxRhs = Float64[]
        auxBeg = Cint[]
        auxNonzIdx = Cint[]
        auxNonzVals = Float64[]
        auxOffset = Cint(0)

        for (k, i) in enumerate(gIntIdx)
            auxCol = Cint(ncols + k - 1)  # 0-based: aux columns follow the original ncols
            xCol = Cint(i - 1)

            # Add the first row for the absolute value constraint: aux_k - x_i >= 0
            push!(auxLhs, 0.0)
            push!(auxRhs, inf)
            push!(auxBeg, auxOffset)
            push!(auxNonzIdx, auxCol)
            push!(auxNonzVals, 1.0)
            push!(auxNonzIdx, xCol)
            push!(auxNonzVals, -1.0)
            auxOffset += 2

            # Add the second row for the absolute value constraint: aux_k + x_i >= 0
            push!(auxLhs, 0.0)
            push!(auxRhs, inf)
            push!(auxBeg, auxOffset)
            push!(auxNonzIdx, auxCol)
            push!(auxNonzVals, 1.0)
            push!(auxNonzIdx, xCol)
            push!(auxNonzVals, 1.0)
            auxOffset += 2
        end

        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(
            lpi,
            Cint(2 * nGInt),
            auxLhs,
            auxRhs,
            C_NULL,
            auxOffset,
            auxBeg,
            auxNonzIdx,
            auxNonzVals
        )
    else
        nGInt = 0
    end

    nlmoCols = ncols + nGInt
    nlmoRows = nrows + 2 * nGInt

    ncolsRef = Ref{Cint}(0)
    nrowsRef = Ref{Cint}(0)
    SCIP.SCIPlpiGetNCols(lpi, ncolsRef)
    SCIP.SCIPlpiGetNRows(lpi, nrowsRef)
    @assert ncolsRef[] == nlmoCols "LPI ncols $(ncolsRef[]) != LMO ncols $nlmoCols"
    @assert nrowsRef[] == nlmoRows "LPI nrows $(nrowsRef[]) != LMO nrows $nlmoRows"

    if verbose > 1
        printstyled("[LMO] ", color=:magenta)
        @printf("Stage %d: ncols=%d (orig=%d aux=%d) nrows=%d (orig=%d aux=%d)\n",
            stage, nlmoCols, ncols, nGInt, nlmoRows, nrows, 2 * nGInt)
    end

    return LPILMO(lpi, nlmoCols, nlmoRows, nrows, verbose, Ref(Inf))
end

# Update the rounding values in the LMO's LPI
function LPIupdateRounding!(lmo::LPILMO, gIntIdx::Vector{Int}, xRound::Vector{Float64})
    nAuxRows = 2 * length(gIntIdx)
    ind = Cint.(lmo.origNrows:(lmo.origNrows + nAuxRows - 1))
    lhs = Float64[]
    rhs = Float64[]
    inf = SCIP.SCIPlpiInfinity(lmo.lpi)

    for i in gIntIdx
        push!(lhs, -xRound[i])
        push!(rhs, inf)
        push!(lhs, xRound[i])
        push!(rhs, inf)
    end

    SCIP.@SCIP_CALL SCIP.SCIPlpiChgSides(lmo.lpi, Cint(nAuxRows), ind, lhs, rhs)
end

# Extract LP basis from SCIP's internal LP solver
function LPIgetBase(scip::Ptr{SCIP.SCIP_}, ncols::Int32, nrows::Int32)
    # Get SCIP's internal LP handle
    lpiPtr = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPgetLPI(scip, lpiPtr)
    lpi = lpiPtr[]

    isOptimal = SCIP.SCIPlpiIsOptimal(lpi)
    @assert isOptimal == SCIP.TRUE "LPI is not optimal"

    # Extract current LP basis
    cstat = zeros(Cint, ncols)
    rstat = zeros(Cint, nrows)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetBase(lpi, cstat, rstat)

    return cstat, rstat
end

# Set LP basis
function LPIsetBase(lmo::LPILMO, cstat::Vector{Cint}, rstat::Vector{Cint})
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetBase(lmo.lpi, cstat, rstat)
end

# Initialize LMO basis from SCIP's internal LP basis
function LPIinitBase(scip::Ptr{SCIP.SCIP_}, lmo::LPILMO, ncols::Int32, nrows::Int32)
    cstatOrig, rstatOrig = LPIgetBase(scip, ncols, nrows)
    if lmo.ncols == ncols && lmo.nrows == nrows
        LPIsetBase(lmo, cstatOrig, rstatOrig)
    else
        cstat = zeros(Cint, lmo.ncols)
        cstat[1:ncols] .= cstatOrig
        rstat = ones(Cint, lmo.nrows)
        rstat[1:nrows] .= rstatOrig
        LPIsetBase(lmo, cstat, rstat)
    end
end

# Compute the extreme point of the LMO's polytope in the direction of the given gradient
# Called in each Frank-Wolfe iteration
function FrankWolfe.compute_extreme_point(lmo::LPILMO, direction::AbstractVector; kwargs...)
    # Set objective
    SCIP.@SCIP_CALL SCIP.SCIPlpiChgObj(
        lmo.lpi,
        Cint(lmo.ncols),
        [Cint(i) for i in 0:(lmo.ncols - 1)],
        Vector{Cdouble}(direction)
    )

    # Capture the basis before solving
    nrowsRef = Ref{Cint}(0)
    SCIP.SCIPlpiGetNRows(lmo.lpi, nrowsRef)
    cstatBefore = zeros(Cint, lmo.ncols)
    rstatBefore = zeros(Cint, nrowsRef[])
    SCIP.SCIPlpiGetBase(lmo.lpi, cstatBefore, rstatBefore)

    # Cap this single solve to whatever's left of the heuristic's own time budget, so one slow LP solve can't run 
    # unbounded and get the whole process killed by an external wall-clock limit instead of exiting cleanly
    lpTimeLeft = min(DEF_SCIP_TIME_LIMIT, max(0.0, lmo.deadline[] - time()))
    SCIP.@SCIP_CALL SCIP.SCIPlpiSetRealpar(lmo.lpi, SCIP.SCIP_LPPAR_LPTILIM, lpTimeLeft)  # SCIP_LPPAR_LPTILIM: LP time limit (> 0)

    # Solve with dual simplex
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolveDual(lmo.lpi)

    # Check if time limit was reached
    if SCIP.SCIPlpiIsTimelimExc(lmo.lpi) == SCIP.TRUE
        throw(LMODeadlineExceeded())
    end

    simplexIter = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetIterations(lmo.lpi, simplexIter)

    # Per-solve basis diagnostics are fine-grained (level 3); skip if no simplex iterations were performed
    verbosed = lmo.verbose >= 2 && simplexIter[] > 0

    if verbosed
        # @printf("LMO direction: pos=%d neg=%d zeros=%d\n",
        #     count(>(0.0), direction), count(<(0.0), direction), count(==(0.0), direction))
        @printf("LMO basis before solve: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d\n",
            count(==(0), cstatBefore), count(==(1), cstatBefore), count(==(2), cstatBefore),
            count(==(0), rstatBefore), count(==(1), rstatBefore), count(==(2), rstatBefore))

        cstat = zeros(Cint, lmo.ncols)
        rstat = zeros(Cint, nrowsRef[])
        SCIP.SCIPlpiGetBase(lmo.lpi, cstat, rstat)

        @printf("LMO basis after solve: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d iters=%d\n",
            count(==(0), cstat), count(==(1), cstat), count(==(2), cstat),
            count(==(0), rstat), count(==(1), rstat), count(==(2), rstat),
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

        if verbosed
            @printf("LMO solution: obj=%.4f\n", obj[])
        end

        return solVector
    else
        error("LPILMO: LP not primal feasible after dual simplex solve")
    end
end