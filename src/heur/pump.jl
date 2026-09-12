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

function roundSolution!(
    xRound::Vector{Float64},
    x::Vector{Float64},
    intIdx::Vector{Int},
    randRound::Bool
)
    xRound .= x  # refresh continuous/inactive-stage entries from the current LP solution
    threshold = getRoundingThreshold(randRound)
    for i in intIdx
        xRound[i] = floor(x[i] + threshold)
    end
end

# Hash function for cycle detection (only hashes integer variable values)
function hashRounded(x::Vector{Float64}, intIdx::Vector{Int})
    hash(tuple((x[i] for i in intIdx)...))
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
function flipVariable!(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
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
    elseif SCIP.SCIPisFeasLT(scip, xRound[i], x[i]) == SCIP.TRUE
        xRound[i] += 1.0
    # Case 4: xRound[i] is an interior point and higher than its fractional point
    elseif SCIP.SCIPisFeasGT(scip, xRound[i], x[i]) == SCIP.TRUE
        xRound[i] -= 1.0
    # Case 5: xRound[i] is an interior point and exactly at its fractional point
    else
        xRound[i] += rand() < 0.5 ? 1.0 : -1.0
    end
end

function perturb(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
    lp::LPInfo,
    st::StageState,
    config::FPFWConfig,
)::Tuple{Int, Bool}
    activeIntIdx = st.activeIntIdx
    nFlips = trunc(Int, st.avgFlips * (rand() + 0.5))
    walksatFired = false

    # Identify the fractional variables in the current rounded solution
    fracVars = Tuple{Float64, Int}[]
    for i in activeIntIdx
        if SCIP.SCIPisFeasEQ(scip, xRound[i], x[i]) == SCIP.FALSE
            roundDist = round(abs(xRound[i] - x[i]), digits=6)
            push!(fracVars, (roundDist, i))
        end
    end

    if isempty(fracVars) && !config.walksatPerturb
        return 0, walksatFired
    end

    # Sort the fractional variables by their distance from the LP solution, descending
    sort!(fracVars, alg=MergeSort, by=first, rev=true)

    if isempty(fracVars)
        nFracFlips = 0
    else
        nFracFlips = clamp(nFlips, 1, length(fracVars))
    end

    selected = fracVars[1:nFracFlips]

    # WalkSAT perturbation: if too few fractional variables to reach nFlips,
    # fill selected with the support of the LP rows that xRound currently violates
    if config.walksatPerturb && nFracFlips < nFlips
        walksatFired = true
        nNeeded = nFlips - nFracFlips
        eligible = st.activeIntIdxSet  # cached, run/stage-invariant - never mutated by infeasibleSupport

        # Get the support of infeasible constraints
        supp = infeasibleSupport(scip, lp.lpRows, xRound, lp.colDict, eligible)
        for (_, i) in selected
            delete!(supp, i)  # delete if already in the selected vars
        end

        xsupp = collect(supp)  # convert to array for shuffling
        # TODO: consider weighting the support by the degree of violation of each constraint,
        # so that more violated constraints are more likely to be selected
        shuffle!(MersenneTwister(config.seed), xsupp)  # MersenneTwister for reproducibility
        for i in xsupp[1:min(nNeeded, length(xsupp))]
            push!(selected, (0.0, i))
        end
    end

    if isempty(selected)
        return 0, walksatFired
    end

    if config.verbose >= 2
        nWalksat = length(selected) - nFracFlips
        println("Perturbing $(length(selected)) integer variables: $nFracFlips fractional, $nWalksat WalkSAT (out of $(length(fracVars)) fractional variables available)")
    end

    # Flip the selected variables
    for (_, i) in selected
        flipVariable!(scip, xRound, x, i, lp.binSet, lp.lpCols)
    end

    return length(selected), walksatFired
end

# Randomize a single general integer within its domain, biased away from its current value once near a bound
function randomizeGeneralInt!(xRound::Vector{Float64}, i::Int, lpCols::Vector{Ptr{SCIP.SCIP_COL}})
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
    x::Vector{Float64},
    prevRound::Vector{Float64},
    lp::LPInfo,
    st::StageState,
    config::FPFWConfig,
    data::FPFWRunData,
)
    activeGIntIdx = st.activeGIntIdx
    avgFlips = st.avgFlips
    changed = 0

    # Binary variables
    for i in lp.binIdx
        r = rand() - 0.47  # [-0.47, 0.53)
        if r > 0 && SCIP.SCIPisFeasEQ(scip, xRound[i], prevRound[i]) == SCIP.TRUE  # stuck variable
            sigma = abs(xRound[i] - x[i])
            if sigma + r > 0.5
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