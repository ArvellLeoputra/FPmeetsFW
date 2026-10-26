# Build the objective function, gradient, and distance function for the Frank-Wolfe projection
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

        f = xExt -> begin
            # |x_j - xRound_j| = (1 - 2 * xRound[j]) * x_j + xRound[j] for binaries x_j
            binDist = sum((1.0 - 2.0 * xRound[i]) * xExt[i] + xRound[i] for i in binIdx; init = 0.0)
            gIntDist = sum(xExt[ncols + k] for k in 1:nGInt; init = 0.0)
            return binDist + gIntDist
        end

        grad! = (storage, xExt) -> begin
            storage[1:ncols] .= 0.0
            for i in binIdx
                storage[i] = 1.0 - 2.0 * xRound[i]
            end
            
            for k in 1:nGInt
                storage[ncols + k] = 1.0
            end
            return storage
        end

        dist = (x, y) -> sum(abs(x[i] - y[i]) for i in intIdx; init = 0.0)

    elseif norm == :euclidean
        f = x -> 0.5 * sum((x[i] - xRound[i])^2 for i in intIdx; init = 0.0)

        grad! = (storage, x) -> begin
            storage .= 0.0
            for i in intIdx
                storage[i] = x[i] - xRound[i]
            end
            return storage
        end

        dist = (x, y) -> 0.5 * sum((x[i] - y[i])^2 for i in intIdx; init = 0.0)
    
    else
        error("Unknown norm: $norm. Choose from :euclidean, :manhattan")
    end

    return f, grad!, dist
end

# Debug-only wrapper around grad!: reports whenever a variable's gradient sign flips between consecutive iterations
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

# Wrap the pure-distance with the objective-FP term:
#   f_alpha(x) = (1-alpha)·f(x) + alpha·weight·sum(objCoeffs[j]·x[j])
# Aux entries for manhattan get the (1-alpha) scale but no objCoef.
function buildOFPFunctions(
    f,
    grad!,
    alpha::Float64,
    weight::Float64,
    objCoeffs::Vector{Float64},
    ncols::Int32
)
    fOFP = x -> begin
        objTerm = 0.0
        for j in 1:ncols
            objTerm += objCoeffs[j] * x[j]
        end
        return (1 - alpha) * f(x) + alpha * weight * objTerm
    end

    gradOFP! = (storage, x) -> begin
        grad!(storage, x)
        storage .*= (1 - alpha)
        for j in 1:ncols
            storage[j] += alpha * weight * objCoeffs[j]
        end
        return storage
    end

    return fOFP, gradOFP!
end

function buildLineSearch(fwStepSize::Symbol, fixedStepSize::Float64=1.0)
    if fwStepSize == :unitary
        FrankWolfe.FixedStep(1.0)
    elseif fwStepSize == :fixed
        FrankWolfe.FixedStep(fixedStepSize)
    elseif fwStepSize == :agnostic
        FrankWolfe.Agnostic()
    elseif fwStepSize == :backtracking
        FrankWolfe.Backtracking()
    elseif fwStepSize == :secant
        FrankWolfe.Secant()
    elseif fwStepSize == :adaptive
        FrankWolfe.Adaptive()
    else
        error("Unknown line search: $fwStepSize. Choose from :unitary, :fixed, :agnostic, :backtracking, :secant, :adaptive")
    end
end

# Callback for runFW: stop a projection early once dual_gap hasn't improved by `minImprovement`
function buildStallCallback(patience::Int, minImprovement::Float64; debug::Bool=false)
    stallCount = Ref(0)
    bestGap = Ref(Inf)
    return (state, args...) -> begin
        gap = state.dual_gap
        if gap < bestGap[] * (1 - minImprovement)
            bestGap[] = gap
            stallCount[] = 0
        else
            stallCount[] += 1
        end
        if debug
            println("    [stall] t=$(state.t) step=$(state.step_type) dual_gap=$gap bestGap=$(bestGap[]) stallCount=$(stallCount[])")
        end
        return stallCount[] < patience
    end
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