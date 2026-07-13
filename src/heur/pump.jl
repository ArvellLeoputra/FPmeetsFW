# Rounding threshold generator
function getRoundingThreshold(random::Bool)
    return random ? rand() : 0.5
end

function roundSolution!(
    xRound::Vector{Float64},
    x::Vector{Float64},
    intIdx::Vector{Int},
    randRound::Bool
)
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
    fracVars = [(abs(xRound[i] - x[i]), i) for i in intIdx if SCIP.SCIPisEQ(scip, xRound[i], x[i]) == SCIP.FALSE]
    if isempty(fracVars)
        # xRound already equals x on every integer variable: there is nothing to perturb.
        return false
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

    return true
end

function restart(
    scip::Ptr{SCIP.SCIP_},
    xRound::Vector{Float64},
    x::Vector{Float64},
    prevRound::Vector{Float64},
    binIdx::Vector{Int},
    gIntIdx::Vector{Int},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    avgFlips::Int
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

    if changed == 0
        for i in binIdx
            if rand() > 0.5
                xRound[i] = 1.0 - xRound[i]
                changed += 1
            end
        end
    end

    # General integer variables
    if !isempty(gIntIdx)
        for _ in 1:avgFlips
            k = rand(1:length(gIntIdx))
            i = gIntIdx[k]

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
            changed += 1
        end
    end
end