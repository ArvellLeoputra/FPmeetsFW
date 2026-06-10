# Function to dispatch Frank-Wolfe variant (for cold-starting)
function callFWVariant(variant::Symbol, f, grad!, lmo, x0; kwargs...)
    if variant == :vanilla
        FrankWolfe.frank_wolfe(f, grad!, lmo, x0; kwargs...)
    elseif variant == :away
        FrankWolfe.away_frank_wolfe(f, grad!, lmo, x0; kwargs...)
    elseif variant == :blended_pairwise
        FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, x0; kwargs...)
    elseif variant == :blended
        FrankWolfe.blended_conditional_gradient(f, grad!, lmo, x0; kwargs...)
    else
        error("Unknown FW variant: $variant. Choose from :vanilla, :away, :blended_pairwise, :blended")
    end
end

# Function to dispatch Frank-Wolfe variant with active set (for warm-starting)
function callFWVariant(variant::Symbol, f, grad!, lmo, activeSet::FrankWolfe.ActiveSet; kwargs...)
    if variant == :vanilla
        error("Vanilla FW variant does not have an active set.")
    elseif variant == :away
        FrankWolfe.away_frank_wolfe(f, grad!, lmo, activeSet; kwargs...)
    elseif variant == :blended_pairwise
        FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, activeSet; kwargs...)
    elseif variant == :blended
        FrankWolfe.blended_conditional_gradient(f, grad!, lmo, activeSet; kwargs...)
    else
        error("Unknown FW variant: $variant. Choose from :away, :blended_pairwise, :blended")
    end
end

function buildFWFunctions(norm::Symbol, intIdx::Vector{Int})
    if norm == :manhattan
        f = (x, xRound) -> sum(abs(x[i] - xRound[i]) for i in intIdx)

        grad! = (storage, x, xRound) -> begin
            storage .= 0.0
            for i in intIdx
                d = x[i] - xRound[i]
                storage[i] = d > 0 ? 1.0 : d < 0 ? -1.0 : 0.0
            end
            return storage
        end

    # check
    elseif norm == :smoothManhattan
        f = (x, xRound) -> sum(sqrt((x[i] - xRound[i])^2 + DEF_INT_TOLERANCE) for i in intIdx)

        grad! = (storage, x, xRound) -> begin
            storage .= 0.0
            for i in intIdx
                d = x[i] - xRound[i]
                storage[i] = d / sqrt(d^2 + DEF_INT_TOLERANCE)
            end
            return storage
        end

    elseif norm == :euclidean
        f = (x, xRound) -> 0.5 * sum((x[i] - xRound[i])^2 for i in intIdx)

        grad! = (storage, x, xRound) -> begin
            storage .= 0.0
            for i in intIdx
                storage[i] = x[i] - xRound[i]
            end
            return storage
        end

    else
        error("Unknown norm: $norm. Choose from :euclidean, :manhattan, :smoothManhattan")
    end

    return f, grad!
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

function run_fw(variant, f, grad!, lmo, x0, activeSet, warmStart, ls, callback, remainingTime)
    if warmStart && activeSet !== nothing && variant !== :vanilla
        callFWVariant(variant, f, grad!, lmo, activeSet,
            max_iteration = DEF_FW_MAX_ITER - 1,
            verbose = false,
            line_search = ls,
            epsilon = DEF_FW_TOLERANCE,
            callback = callback,
            timeout = remainingTime
        )
    else
        callFWVariant(variant, f, grad!, lmo, x0,
            max_iteration = DEF_FW_MAX_ITER - 1,
            verbose = false,
            line_search = ls,
            epsilon = DEF_FW_TOLERANCE,
            callback = callback,
            timeout = remainingTime
        )
    end
end