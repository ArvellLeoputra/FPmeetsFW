function buildGradFunction(intIdx::Vector{Int}, gradFormula)
    return (storage, x) -> begin
        storage .= 0.0
        for i in intIdx
            storage[i] = gradFormula(i, x)
        end
        return storage
    end
end

function buildFWFunctions(norm::Symbol, intIdx::Vector{Int}, xRound::Vector{Float64})
    if norm == :manhattan
        f = x -> sum(abs(x[i] - xRound[i]) for i in intIdx)

        # Currently, only works for binary
        # TODO: General integer case
        grad! = buildGradFunction(intIdx, (i, x) -> begin
            d = x[i] - xRound[i]
            if d > DEF_INT_TOLERANCE
                1.0
            elseif d < -DEF_INT_TOLERANCE
                -1.0
            else
                xRound[i] < 0.5 ? 1.0 : -1.0
            end
        end)

        dist = (x, y) -> sum(abs(x[i] - y[i]) for i in intIdx)
    elseif norm == :smoothManhattan
        f = x -> sum(sqrt((x[i] - xRound[i])^2 + DEF_INT_TOLERANCE) for i in intIdx)

        grad! = buildGradFunction(intIdx, (i, x) -> begin
            d = x[i] - xRound[i]
            d / sqrt(d^2 + DEF_INT_TOLERANCE)
        end)

        dist = (x, y) -> sum(sqrt((x[i] - y[i])^2 + DEF_INT_TOLERANCE) for i in intIdx)
    elseif norm == :euclidean
        f = x -> 0.5 * sum((x[i] - xRound[i])^2 for i in intIdx)

        grad! = buildGradFunction(intIdx, (i, x) -> x[i] - xRound[i])

        dist = (x, y) -> 0.5 * sum((x[i] - y[i])^2 for i in intIdx)
    else
        error("Unknown norm: $norm. Choose from :euclidean, :manhattan, :smoothManhattan")
    end

    return f, grad!, dist
end

function buildLineSearch(lineSearch::Symbol)
    if lineSearch == :unitary
        FrankWolfe.FixedStep(1.0)
    elseif lineSearch == :agnostic
        FrankWolfe.Agnostic()
    elseif lineSearch == :backtracking
        FrankWolfe.Backtracking()
    elseif lineSearch == :secant
        FrankWolfe.Secant()
    elseif lineSearch == :adaptive
        FrankWolfe.Adaptive()
    else
        error("Unknown line search: $lineSearch. Choose from :agnostic, :backtracking, :secant, :adaptive")
    end
end

function runFW(variant, f, grad!, lmo; x0, activeSet, warmStart, ls, remainingTime, fwMaxIterations, callback=nothing, verbose=false)
    if warmStart && activeSet !== nothing && variant !== :vanilla
        startPoint = activeSet
    else
        startPoint = x0
    end

    callFWVariant(variant, f, grad!, lmo, startPoint,
        max_iteration = fwMaxIterations - 1,
        verbose = verbose,
        line_search = ls,
        epsilon = DEF_FW_TOLERANCE,
        callback = callback,
        timeout = remainingTime,
        dual_gap_compute_frequency = 1,
        trajectory = true
    )
end

# Dispatch to the right FrankWolfe.jl algorithm. startPoint is either a plain
# vector (cold-starting) or a FrankWolfe.ActiveSet (warm-starting).
function callFWVariant(variant::Symbol, f, grad!, lmo, startPoint; kwargs...)
    if variant == :vanilla
        if startPoint isa FrankWolfe.ActiveSet
            error("Vanilla FW variant does not have an active set.")
        end
        FrankWolfe.frank_wolfe(f, grad!, lmo, startPoint; kwargs...)
    elseif variant == :away
        FrankWolfe.away_frank_wolfe(f, grad!, lmo, startPoint; kwargs...)
    elseif variant == :blended_pairwise
        FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, startPoint; kwargs...)
    elseif variant == :blended
        FrankWolfe.blended_conditional_gradient(f, grad!, lmo, startPoint; kwargs...)
    else
        error("Unknown FW variant: $variant. Choose from :vanilla, :away, :blended_pairwise, :blended")
    end
end