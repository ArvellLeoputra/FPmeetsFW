# Rounding threshold generator (from Bertacco et al., 2007)
function getRoundingThreshold(random::Bool)
    if !random
        return 0.5
    end

    omega = rand()
    if omega <= 0.5
        return 2 * omega * (1 - omega)
    else
        return 1 - 2 * omega * (1 - omega)
    end
end

# Weighted-random pick of one active-set atom, probability proportional to its FW weight
function sampleActiveSetAtom(activeSet)
    r = rand() * sum(activeSet.weights)
    cum = 0.0
    for i in eachindex(activeSet.weights)
        cum += activeSet.weights[i]
        if r <= cum
            return activeSet.atoms[i]
        end
    end
    return activeSet.atoms[end]  # floating-point fallback
end

function roundSolution!(
    xRound::Vector{Float64},
    x::Vector{Float64},
    intIdx::Vector{Int},
    randRound::Bool
)::Nothing
    xRound .= x  # refresh continuous/inactive-stage entries from the current LP solution
    threshold = getRoundingThreshold(randRound)
    for i in intIdx
        xRound[i] = floor(x[i] + threshold)
    end
end

# Hash function for cycle detection (only hashes integer variable values).
# Hashes the already-integer-valued Float64s directly rather than converting
# to Int, since a diverged iterate can round to a value outside Int64's range.
function hashRounded(xRound::Vector{Float64}, intIdx::Vector{Int})
    return hash(@view xRound[intIdx])
end

# Compute the support of the violated constraints in the LP, restricted to eligible variables
function infeasibleSupport(
    scip::Ptr{SCIP.SCIP_},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    sol::Vector{Float64},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    eligible::Set{Int},
)::Set{Int}
    supp = Set{Int}()
    inf = SCIP.SCIPinfinity(scip)

    for row in lpRows
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

        # Check if the constraint is violated
        belowLhs = lhs > -inf && SCIP.SCIPisFeasLT(scip, activity, lhs) == SCIP.TRUE
        aboveRhs = rhs < inf && SCIP.SCIPisFeasGT(scip, activity, rhs) == SCIP.TRUE

        # Find the support of the violated constraint, but only include eligible variables
        if belowLhs || aboveRhs
            for k in 1:nnonz
                j = colDict[nonzCols[k]]
                if j in eligible
                    push!(supp, j)
                end
            end
        end
    end

    return supp
end

# Flip a single perturbation candidate
function flipRounded!(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    xFrac::Vector{Float64},
    i::Int,
    binSet::Set{Int},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}}
)
    if i in binSet
        xRound[i] = 1.0 - xRound[i]
        return
    end

    var = SCIP.SCIPcolGetVar(lpCols[i])
    lb = SCIP.SCIPvarGetLbLocal(var)
    ub = SCIP.SCIPvarGetUbLocal(var)

    # Case 1: xRound[i] is at the lower bound
    if SCIP.SCIPisFeasEQ(scip, xRound[i], lb) == SCIP.TRUE
        xRound[i] += 1.0
    # Case 2: xRound[i] is at the upper bound
    elseif SCIP.SCIPisFeasEQ(scip, xRound[i], ub) == SCIP.TRUE
        xRound[i] -= 1.0
    # Case 3: xRound[i] is an interior point and lower than its fractional point
    elseif SCIP.SCIPisFeasLT(scip, xRound[i], xFrac[i]) == SCIP.TRUE
        xRound[i] += 1.0
    # Case 4: xRound[i] is an interior point and higher than its fractional point
    elseif SCIP.SCIPisFeasGT(scip, xRound[i], xFrac[i]) == SCIP.TRUE
        xRound[i] -= 1.0
    # Case 5: xRound[i] is an interior point and exactly at its fractional point
    else
        xRound[i] += rand() < 0.5 ? 1.0 : -1.0
    end
end

function perturb(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    xFrac::Vector{Float64},
    lp::LPInfo,
    st::StageState,
    config::FPFWConfig,
)::Tuple{Int, Int}
    activeIntIdx = st.activeIntIdx
    nFlips = trunc(Int, st.avgFlips * (rand() + 0.5))

    # Identify the fractional variables in the current rounded solution
    fracVars = Tuple{Float64, Int}[]
    for i in activeIntIdx
        if SCIP.SCIPisFeasEQ(scip, xRound[i], xFrac[i]) == SCIP.FALSE
            roundDist = round(abs(xRound[i] - xFrac[i]), digits=6)
            push!(fracVars, (roundDist, i))
        end
    end

    if isempty(fracVars) && !config.walksatPerturb
        return 0, 0
    end

    # Sort the fractional variables by their distance from the LP solution, descending
    # MergeSort is used due to stability, so equal distances keep ascending-index order
    sort!(fracVars, alg=MergeSort, by=first, rev=true)

    # Calculate number of fractional variables to flip
    if isempty(fracVars)
        nFracFlips = 0
    else
        nFracFlips = clamp(nFlips, 1, length(fracVars))
    end

    # WalkSAT perturbation: if there are too few fractional variables to reach the flip target,
    # fill the rest with a random subset of the variables appearing in violated constraints.
    varsToFlip = fracVars[1:nFracFlips]
    if config.walksatPerturb && nFracFlips < nFlips
        nNeeded = nFlips - nFracFlips
        eligible = st.activeIntIdxSet  # cached, run/stage-invariant - never mutated by infeasibleSupport

        # Get the support of infeasible constraints
        supp = infeasibleSupport(scip, lp.lpRows, xRound, lp.colDict, eligible)
        for (_, i) in fracVars
            delete!(supp, i)  # delete if already in the selected vars
        end

        xsupp = collect(supp)  # convert to array for shuffling
        # TODO:  weight shuffle the variables, the more the variables appear inside violated constraints, the "hopefully" better it is to perturb this
        shuffle!(MersenneTwister(config.seed), xsupp)  # MersenneTwister for reproducibility
        for i in xsupp[1:min(nNeeded, length(xsupp))]
            push!(varsToFlip, (0.0, i))
        end
    end

    if isempty(varsToFlip)
        return 0, 0
    end

    nWalksat = length(varsToFlip) - nFracFlips

    if config.verbose >= 2
        println("Perturbing $(length(varsToFlip)) integer variables: $nFracFlips fractional, $nWalksat WalkSAT (out of $(length(fracVars)) fractional variables available)")
    end

    # Flip the selected variables
    for (_, i) in varsToFlip
        flipRounded!(scip, xRound, xFrac, i, lp.binSet, lp.lpCols)
    end

    return length(varsToFlip), nWalksat
end

# Randomize a single general integer within its domain, biased away from its current value once near a bound
function randomizeGeneralInt!(
    xRound::Vector{Float64},
    i::Int,
    lpCols::Vector{Ptr{SCIP.SCIP_COL}}
)
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
        xRound[i] + (2 * DEF_BIGM - 1) * r - DEF_BIGM
    end

    xRound[i] = clamp(floor(newVal), lb, ub)
end

function restart(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    xFrac::Vector{Float64},
    prevRound::Vector{Float64},
    lp::LPInfo,
    st::StageState,
    config::FPFWConfig,
    data::FPFWRunData
)::Int
    activeGIntIdx = st.activeGIntIdx
    avgFlips = st.avgFlips
    changed = 0

    # Binary variables
    for i in lp.binIdx
        r = rand() - 0.47  # [-0.47, 0.53)
        if r > 0 && SCIP.SCIPisFeasEQ(scip, xRound[i], prevRound[i]) == SCIP.TRUE  # stuck variable
            roundDist = abs(xRound[i] - xFrac[i])
            if roundDist + r > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    # General integer variables (stage 2 only)
    if !isempty(activeGIntIdx)
        # Count the number of iterations since the last restart
        nitr = data.stats.pumpIterations
        staleIters = nitr - data.lastRestart

        # Apply geometric decay to flip budget
        for _ in 1:staleIters
            data.gIntFlipBudget = trunc(Int, data.gIntFlipBudget * DEF_GEOM_FACTOR)
        end
        data.lastRestart = nitr

        # Increase the flip budget by 2 * avgFlips + 1, but cap it at gIntFlipBudgetCap
        data.gIntFlipBudget = min(data.gIntFlipBudget + 2 * avgFlips + 1, data.gIntFlipBudgetCap)

        for _ in 1:data.gIntFlipBudget
            # Randomly select a general integer variable to randomize
            i = activeGIntIdx[rand(1:length(activeGIntIdx))]
            prevVal = xRound[i]
            randomizeGeneralInt!(xRound, i, lp.lpCols)

            # Check if the variable's value has changed
            if SCIP.SCIPisFeasEQ(scip, xRound[i], prevVal) == SCIP.FALSE
                changed += 1
            end
        end
    else
        # Safety net: no general integers to perturb, so shake all binaries hard if nothing changed above
        if changed == 0
            for i in lp.binIdx
                if rand() > 0.5
                    xRound[i] = 1.0 - xRound[i]
                    changed += 1
                end
            end
        end
    end

    if config.verbose >= 2
        println("Restarting: $changed variables changed")
    end

    return changed
end

# Randomized feasibility probe: probabilistically round the current LP point xFrac several times and
# submit each to SCIP. Returns true as soon as one is accepted, leaving the MIP-feasible point in xProbe.
# Called at both stages with all integer constrained variables considered.
function randFeasCheck!(
    scip::Ptr{SCIP.SCIP_},
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
    lp::LPInfo,
    xFrac::Vector{Float64},
    xProbe::Vector{Float64}
)::Bool
    k = count(i -> !isVarInteger(scip, xFrac[i]), lp.intIdx)

    # Cap the nFrac to minimize randRound time
    if k == 0 || k > DEF_RAND_FEAS_MAX_FRAC
        return false
    end

    # Cap nProbes to 2^k, where k = # fractional integer vars
    nProbes = min(DEF_RAND_FEAS_ITER_LIMIT, 1 << min(k, 10))

    for _ in 1:nProbes
        xProbe .= xFrac
        for i in lp.intIdx
            # Round up with probability xFrac[i] - floor(xFrac[i]), otherwise round down
            frac = xFrac[i] - floor(xFrac[i])
            xProbe[i] = rand() < frac ? ceil(xFrac[i]) : floor(xFrac[i])
        end

        if submitSolution(scip, heur_ptr, lp.lpCols, xProbe, lp.ncols)
            return true
        end
    end
    return false
end