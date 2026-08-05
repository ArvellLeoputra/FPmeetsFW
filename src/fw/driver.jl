# Build a gradient function for the specified indices and formula
function buildGradFunction(intIdx::Vector{Int}, gradFormula)
    return (storage, x) -> begin
        storage .= 0.0
        for i in intIdx
            storage[i] = gradFormula(i, x)
        end
        return storage
    end
end

# Debug-only wrapper around grad!: reports whenever a variable's gradient sign flips between
# consecutive evaluations (prevGrad holds the previous gradient). Only used at verbose >= 2.
function buildGradCheck(grad!, prevGrad, intIdx, xRound)
    return (storage, x) -> begin
        grad!(storage, x)
        for i in intIdx
            if abs(x[i] - xRound[i]) > DEF_INT_TOLERANCE &&
                abs(prevGrad[i]) > DEF_INT_TOLERANCE && 
                sign(storage[i]) != sign(prevGrad[i])
                @printf("  [grad flip] var=%d x=%.6f xRound=%.6f old=%.4f new=%.4f\n",
                    i, x[i], xRound[i], prevGrad[i], storage[i])
            end
        end
        copyto!(prevGrad, storage)
        return storage
    end
end

function buildFWFunctions(
    norm::Symbol,
    binIdx::Vector{Int},
    gIntIdx::Vector{Int},
    intIdx::Vector{Int},
    xRound::Vector{Float64}
)
    if norm == :manhattan
        ncols = length(xRound)
        nGInt = length(gIntIdx)

        function binTerm(i, x)
            if xRound[i] < 0.5
                return x[i]
            else
                return 1.0 - x[i]
            end
        end

        f = xExt -> begin
            binDist = sum(binTerm(i, xExt) for i in binIdx; init = 0.0)
            gIntDist = sum(xExt[ncols + k] for k in 1:nGInt; init = 0.0)
            binDist + gIntDist
        end

        grad! = (storage, xExt) -> begin
            storage .= 0.0
            for i in binIdx
                storage[i] = xRound[i] < 0.5 ? 1.0 : -1.0
            end

            for k in 1:nGInt
                storage[ncols + k] = 1.0
            end
            return storage
        end

        dist = (x, y) -> sum(abs(x[i] - y[i]) for i in intIdx; init = 0.0)

    elseif norm == :smoothManhattan
        f = x -> sum(sqrt((x[i] - xRound[i])^2 + DEF_INT_TOLERANCE) for i in intIdx; init = 0.0)

        grad! = buildGradFunction(intIdx, (i, x) -> begin
            d = x[i] - xRound[i]
            d / sqrt(d^2 + DEF_INT_TOLERANCE)
        end)

        dist = (x, y) -> sum(sqrt((x[i] - y[i])^2 + DEF_INT_TOLERANCE) for i in intIdx; init = 0.0)
    
    elseif norm == :euclidean
        f = x -> 0.5 * sum((x[i] - xRound[i])^2 for i in intIdx; init = 0.0)

        grad! = buildGradFunction(intIdx, (i, x) -> x[i] - xRound[i])

        dist = (x, y) -> 0.5 * sum((x[i] - y[i])^2 for i in intIdx; init = 0.0)
    
    else
        error("Unknown norm: $norm. Choose from :euclidean, :manhattan, :smoothManhattan")
    end

    return f, grad!, dist
end

function buildLineSearch(fwStepSize::Symbol)
    if fwStepSize == :unitary
        FrankWolfe.FixedStep(1.0)
    elseif fwStepSize == :agnostic
        FrankWolfe.Agnostic()
    elseif fwStepSize == :backtracking
        FrankWolfe.Backtracking()
    elseif fwStepSize == :secant
        FrankWolfe.Secant()
    elseif fwStepSize == :adaptive
        FrankWolfe.Adaptive()
    else
        error("Unknown line search: $fwStepSize. Choose from :unitary, :agnostic, :backtracking, :secant, :adaptive")
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
        trajectory = true
    )
end

function callFWVariant(variant::Symbol, f, grad!, lmo, startPoint; kwargs...)
    if variant == :vanilla
        if startPoint isa FrankWolfe.ActiveSet
            error("Vanilla FW variant does not have an active set.")
        end
        FrankWolfe.frank_wolfe(f, grad!, lmo, startPoint; dual_gap_compute_frequency=1, kwargs...)
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