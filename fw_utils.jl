# Function to dispatch Frank-Wolfe variant (for cold-starting)
function call_fw_variant(variant::Symbol, f, grad!, lmo, x0; kwargs...)
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
function call_fw_variant(variant::Symbol, f, grad!, lmo, active_set::FrankWolfe.ActiveSet; kwargs...)                                                                                                                             
    if variant == :vanilla                                                                                                                                                                                                      
        error("Vanilla FW variant does not have an active set.")
    elseif variant == :away                                                                                                                                                                                                       
        FrankWolfe.away_frank_wolfe(f, grad!, lmo, active_set; kwargs...)
    elseif variant == :blended_pairwise                                                                                                                                                                                           
        FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, active_set; kwargs...)                                                                                                                                  
    elseif variant == :blended                                                                                                                                                                                                    
        FrankWolfe.blended_conditional_gradient(f, grad!, lmo, active_set; kwargs...)
    else                                                                                                                                                                                                                          
        error("Unknown FW variant: $variant. Choose from :away, :blended_pairwise, :blended")                                                                                                                         
    end                                                                                                                                                                                                                           
end

function build_fw_functions(norm::Symbol, all_integers::Vector{Int})
    if norm == :manhattan
        f = (x, x_round) -> sum(abs(x[i] - x_round[i]) for i in all_integers)

        grad! = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in all_integers
                d = x[i] - x_round[i]
                storage[i] = d > 0 ? 1.0 : d < 0 ? -1.0 : 0.0
            end
            return storage
        end

    elseif norm == :smooth_manhattan
        f = (x, x_round) -> sum(sqrt((x[i] - x_round[i])^2 + DEF_TOLERANCE) for i in all_integers)

        grad! = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in all_integers
                d = x[i] - x_round[i]
                storage[i] = d / sqrt(d^2 + DEF_TOLERANCE)
            end
            return storage
        end

    elseif norm == :euclidean
        f = (x, x_round) -> 0.5 * sum((x[i] - x_round[i])^2 for i in all_integers)

        grad! = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in all_integers
                storage[i] = x[i] - x_round[i]
            end
            return storage
        end

    else
        error("Unknown norm: $norm. Choose from :euclidean, :manhattan, :smooth_manhattan")
    end

    return f, grad!
end

function build_line_search(line_search::Symbol)
    if line_search == :agnostic
        FrankWolfe.Agnostic()
    elseif line_search == :backtracking
        FrankWolfe.Backtracking()
    elseif line_search == :secant
        FrankWolfe.Secant()
    elseif line_search == :adaptive
        FrankWolfe.Adaptive()
    else
        error("Unknown line search: $line_search. Choose from :agnostic, :backtracking, :secant, :adaptive")
    end
end

function run_fw(variant, f, grad!, lmo, x0, active_set, warm_start, ls, callback, remaining_time)
    if warm_start && active_set !== nothing && variant !== :vanilla
        call_fw_variant(variant, f, grad!, lmo, active_set,
            max_iteration = DEF_FW_MAX_ITER,
            verbose = false,
            line_search = ls,
            epsilon = DEF_FW_TOLERANCE,
            callback = callback,
            timeout = remaining_time
        )
    else
        call_fw_variant(variant, f, grad!, lmo, x0,
            max_iteration = DEF_FW_MAX_ITER,
            verbose = false,
            line_search = ls,
            epsilon = DEF_FW_TOLERANCE,
            callback = callback,
            timeout = remaining_time
        )
    end
end