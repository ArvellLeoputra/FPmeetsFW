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

function perturb(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
    binIdx::Vector{Int},
    intIdx::Vector{Int},
    avgFlips::Int,
    verbose::Bool
)
    fracVars = [(round(abs(xRound[i] - x[i]), digits=6), i) for i in intIdx if SCIP.SCIPisEQ(scip, xRound[i], x[i]) == SCIP.FALSE]

    if isempty(fracVars)
        return 0
    end
    sort!(fracVars, rev=true)

    nFlips = round(Int, avgFlips * (rand() + 0.5))
    nFracFlips = clamp(nFlips, 1, length(fracVars))

    if verbose
        println("Perturbing $nFracFlips integer variables (out of $(length(fracVars)) fractional variables)")
    end

    binSet = Set(binIdx)  # faster lookup
    for k in 1:nFracFlips
        _, i = fracVars[k]
        if i in binSet  # binary: round to opposite
            xRound[i] = 1.0 - xRound[i]
        else  # general integer: reverse rounding direction
            xRound[i] += SCIP.SCIPisLT(scip, xRound[i], x[i]) == SCIP.TRUE ? 1.0 : -1.0
        end
    end

    return nFracFlips
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
    binIdx::Vector{Int},
    gIntIdx::Vector{Int},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    avgFlips::Int,
    verbose::Bool
)
    changed = 0

    # Binary variables
    for i in binIdx
        r = rand() - 0.47  # [-0.47, 0.53)
        if r > 0 && SCIP.SCIPisEQ(scip, xRound[i], prevRound[i]) == SCIP.TRUE  # stuck variable
            sigma = abs(xRound[i] - x[i])
            if sigma + r > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    # General integer variables (stage 2 only)
    # TODO: compare the fixed avgFlips to an adaptive count that grows with the number of restarts
    if !isempty(gIntIdx)
        for _ in 1:avgFlips
            i = gIntIdx[rand(1:length(gIntIdx))]
            prevVal = xRound[i]
            randomizeGeneralInt!(xRound, i, lpCols)
            if SCIP.SCIPisEQ(scip, xRound[i], prevVal) == SCIP.FALSE
                changed += 1
            end
        end
    end

    # Safety net: only shake binaries aggressively if nothing changed above
    if changed == 0
        for i in binIdx
            if rand() > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    if verbose
        println("Restarting: $changed variables changed")
    end

    return changed
end