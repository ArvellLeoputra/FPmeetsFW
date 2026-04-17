# Function to assign Frank-Wolfe variant.
function call_fw_variant(variant::Symbol, f, grad!, lmo, x0; kwargs...)
    if variant == :vanilla
        return FrankWolfe.frank_wolfe(f, grad!, lmo, x0; kwargs...)
    elseif variant == :away
        return FrankWolfe.away_frank_wolfe(f, grad!, lmo, x0; kwargs...)
    elseif variant == :blended_pairwise
        return FrankWolfe.blended_pairwise_conditional_gradient(f, grad!, lmo, x0; kwargs...)
    elseif variant == :blended
        return FrankWolfe.blended_conditional_gradient(f, grad!, lmo, x0; kwargs...)
    else
        error("Unknown FW variant: $variant. Choose from :vanilla, :away, :blended_pairwise, :blended")
    end
end

# Main FPFW Heuristic Implementation
function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}

    @assert SCIP.SCIPhasCurrentNodeLP(scip) == SCIP.TRUE  # always true, since we set timing to DURINGLPLOOP
    result = SCIP.SCIP_DIDNOTFIND

    heur.called += 1
    if heur.called > 1
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    if DEF_RANDOM_SEED !== nothing
        Random.seed!(DEF_RANDOM_SEED)
    end

    stats = FPFWStats()
    total_start_time = time()

    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)

    # Build LMO from current LP
    if heur.lmo === nothing
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
    end

    lmo = heur.lmo

    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)
    col_to_idx = Dict(lp_cols[k] => k for k in 1:nvars)

    binary = Int[]
    current_solution = zeros(SCIP.SCIP_Real, nvars)

    # Identify binary variables and get current LP solution
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)

        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binary, j)
        end

        current_solution[j] = SCIP.SCIPcolGetPrimsol(col) # Get current LP solution
    end

    if DEBUG_VERBOSE
        println("Initial solution:")
        for j in 1:nvars
            @printf("  x[%d] = %.3f\n", j, current_solution[j])
        end
    end

    lp_obj = sum(current_solution[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
    stats.lp_objective = lp_obj
    
    if DEBUG_VERBOSE
        println("Initial LP objective: $(round(lp_obj, digits=4))")
        n_frac = 0
        for i in binary
            if abs(current_solution[i] - round(current_solution[i])) > DEF_TOLERANCE
                n_frac += 1
            end
        end
        println("Non-integral binary vars: $n_frac/$(length(binary))")
    end

    if heur.projection_norm == :manhattan
        f_proj = (x, x_round) -> begin
            sum = 0.0
            for i in binary
                sum += abs(x[i] - x_round[i])
            end
            return sum
        end

        g_proj = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in binary
                d = x[i] - x_round[i]
                if d > 0
                    storage[i] = 1.0
                elseif d < 0
                    storage[i] = -1.0
                else
                    storage[i] = 0.0
                end
            end
            return storage
        end

    elseif heur.projection_norm == :euclidean
        f_proj = (x, x_round) -> begin
            sum = 0.0
            for i in binary
                d = x[i] - x_round[i]
                sum += d * d
            end
            return 0.5 * sum
        end

        g_proj = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in binary
                storage[i] = x[i] - x_round[i]
            end
            return storage
        end

    elseif heur.projection_norm == :abssmooth
        f_proj = (x, x_round) -> begin
            sum = 0.0
            for i in binary
                d = x[i] - x_round[i]
                sum += sqrt(d * d + DEF_TOLERANCE)  # Smooth approximation of abs
            end
            return sum
        end

        g_proj = (storage, x, x_round) -> begin
            storage .= 0.0
            for i in binary
                d = x[i] - x_round[i]
                storage[i] = d / sqrt(d * d + DEF_TOLERANCE)  # Gradient of smooth abs
            end
            return storage
        end

    else
        error("Unknown projection norm: $(heur.projection_norm)")
    end

    # Hash function for cycle detection (only hashes binary variable values)
    function hash_rounded(x_round, binary_indices)
        return hash(tuple((x_round[i] for i in binary_indices)...))
    end

    # Rounding function using custom threshold
    function round_with_threshold(x, binary_indices, threshold)
        x_round = zeros(Float64, length(x))

        for i in binary_indices
            x_round[i] = x[i] >= threshold ? 1.0 : 0.0
        end

        return x_round
    end

    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(current_solution)  # Initial LP feasible solution
    prev_x = copy(x)  # To calculate distance moved
    last_projected_x = copy(current_solution)  # Always keep a feasible point for restarts
    x_after = zeros(Float64, nvars)  # Preallocated buffer for post-step x in callback

    # Cycle detection: store hashes of all visited rounded solutions
    rounded_cache = Set{UInt64}()
    fw_traj = Vector{Tuple{Int, Vector{Float64}, Float64}}()

    # TODO: store the best solution found across iterations, not just the first one
    found_solution = nothing

    ls = if heur.line_search == :backtracking
        FrankWolfe.Backtracking()
    elseif heur.line_search == :secant
        FrankWolfe.Secant()
    elseif heur.line_search == :adaptive
        FrankWolfe.Adaptive()
    else
        FrankWolfe.Agnostic()
    end

    iter = 0
    while true
        iter += 1
        if DEBUG_VERBOSE
            println("\n--- FPFW Iteration $iter ---")
        end
        stats.fp_iterations = iter

        # Check time limit
        elapsed = time() - total_start_time
        if elapsed > DEF_FW_TIME_LIMIT
            println("FW time limit reached ($(round(elapsed, digits=2))s > $(DEF_FW_TIME_LIMIT)s)")
            break
        end
        remaining_time = DEF_FW_TIME_LIMIT - elapsed

        # Step 1: Rounding LP feasible solution with custom threshold
        x_round = round_with_threshold(x, binary, heur.rounding_threshold)

        if DEBUG_VERBOSE
            println("Rounding(threshold=$(heur.rounding_threshold)):")
            for i in binary
                @printf("  x[%d]: %.3f -> %.1f\n", i, x[i], x_round[i])
            end
            println("------------------------ \n")
        end

        # Cycle detection: check if we've visited this rounded solution before
        h = hash_rounded(x_round, binary)
        if h in rounded_cache
            stats.restart_cycles += 1
            if DEBUG_VERBOSE
                println("Cycle detected at iteration $iter (restart #$(stats.restart_cycles))")
            end

            if stats.restart_cycles >= DEF_MAX_CYCLE_RESTARTS
                println("Maximum cycle restarts ($DEF_MAX_CYCLE_RESTARTS) reached, stopping.")
                break
            end

            # Perturb the rounding target
            if DEBUG_VERBOSE
                println("Perturbing:")
            end

            for i in binary
                old = x_round[i]

                if rand() < DEF_PERTURB_FRACTION
                    x_round[i] = 1.0 - x_round[i]
                end
                
                if DEBUG_VERBOSE
                    @printf("  x[%d]: %.1f -> %.1f\n", i, old, x_round[i])
                end
            end

            # Reset x to last known feasible point for starting the next FW iteration
            x .= last_projected_x
            prev_x .= last_projected_x

            if DEBUG_VERBOSE
                println("restarting from feasible point... \n")
            end

            # Recompute hash after perturbation so we cache the new rounding
            h = hash_rounded(x_round, binary)
        end

        push!(rounded_cache, h)

        empty!(fw_traj)
        fw_callback = (state, args...) -> begin
            if state.step_type !== FrankWolfe.ST_LAST && state.step_type !== FrankWolfe.ST_POSTPROCESS
                x_after .= state.x .- state.gamma .* state.d
                push!(fw_traj, (state.t, copy(x_after), state.primal))
            end
            return true
        end

        # Step 2: Projection using Frank-Wolfe
        # Line search selection
        fw_start = time()
        x_new, _ = call_fw_variant(
            heur.fw_variant,
            x -> f_proj(x, x_round),
            (storage, x) -> g_proj(storage, x, x_round),
            lmo,
            x,
            max_iteration = DEF_FW_MAX_ITER,
            verbose = false,
            line_search = ls,
            epsilon = DEF_FW_TOLERANCE,
            callback = fw_callback,
            timeout = remaining_time
        )
        stats.fw_time += time() - fw_start

        fw_iters = length(fw_traj)
        stats.fw_iterations += fw_iters

        if DEBUG_VERBOSE
            for (t, xk, fk) in fw_traj
                @printf("--- FW Step %d --- \n", t)
                @printf("Solution: \n")
                for i in binary
                    @printf("  x[%d]: %.3f \n", i, xk[i])
                end
                @printf("Objective = %.3f\n \n", fk)
            end
        end

        # Step 3: Check feasibility and integrality and calculate current objective value
        is_integral = check_integrality(x_new, binary)
        is_feasible = check_feasibility(scip, x_new, col_to_idx)
        distance = sqrt(sum((x_new[i] - prev_x[i])^2 for i in binary))
        obj_val = sum(x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)

        if DEBUG_VERBOSE
            n_frac = 0
            for i in binary
                if abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE
                    n_frac += 1
                end
            end

            total_dist = sum(abs(x_new[i] - x_round[i]) for i in binary)
            
            @printf("[Iter %d] obj=%.6f | dist_moved=%.6f | dist_target=%.6f | frac=%d/%d | fw_iters=%d | integral=%s | feasible=%s\n",
                iter, obj_val, distance, total_dist, n_frac, length(binary), fw_iters, is_integral, is_feasible)
        end

        # Just for debugging (numerical errors)
        if !is_feasible
            println("FrankWolfe return infeasible solution!")
            stats.infeasible_exit = true
            break
        end

        if is_integral
            # Submit solution to SCIP
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_new, nvars)
                stats.solution_found = true
                stats.iter_found_solution = iter
                stats.final_objective = obj_val
                result = SCIP.SCIP_FOUNDSOL

                if DEBUG_VERBOSE
                    println("Solution:")
                    for j in 1:nvars
                        @printf("  x[%d] = %.3f\n", j, x_new[j])
                    end
                end
                break
            else
                println("Solution rejected by SCIP")
                break
            end
        else
            # Fixed-point detection: solution barely changed
            if distance < DEF_TOLERANCE
                stats.restart_fixedpoint += 1
                if DEBUG_VERBOSE
                    println("\n Converged to fixed point at iteration $iter (restart #$(stats.restart_fixedpoint))")
                end

                if stats.restart_fixedpoint >= DEF_MAX_FIXEDPOINT_RESTARTS
                    println("Maximum fixed-point restarts ($DEF_MAX_FIXEDPOINT_RESTARTS) reached, stopping.")
                    break
                end

                # Reset to feasible point and let next iteration create a new rounding target
                x .= last_projected_x
                prev_x .= last_projected_x

                # Slightly perturb x within feasible bounds to get different rounding
                if DEBUG_VERBOSE
                    println("Perturbing:")
                end

                for i in binary
                    old = x[i]
                    x[i] = clamp(x[i] + DEF_FIXEDPOINT_PERTURB * (rand() - 0.5), 0.0, 1.0)

                    if DEBUG_VERBOSE
                        @printf("  x[%d]: %.3f -> %.3f\n", i, old, x[i])
                    end
                end

                if DEBUG_VERBOSE
                    println("restarting from feasible point... \n")
                end
                continue
            end

            # Continue with the projected solution for next FW iteration
            prev_x .= x
            x .= x_new
            last_projected_x .= x_new  # Store as last known feasible point
        end

    end

    # Print final summary
    total_time = time() - total_start_time
    stats.total_time = total_time
    println("\n" * "="^80)
    println("FPFW HEURISTIC SUMMARY")
    println("="^80)
    println("Binary variables:  $(length(binary))")
    println("Total time:        $(round(total_time, digits=2)) seconds")
    println("FP iterations:     $(stats.fp_iterations)")
    println("FW iterations:     $(stats.fw_iterations)")
    println("FW time:           $(round(stats.fw_time, digits=2)) seconds")
    println("Restarts (cycles): $(stats.restart_cycles)")
    println("Restarts (fixed):  $(stats.restart_fixedpoint)")
    println("Solution found:    $(stats.solution_found)")

    if stats.solution_found
        println("Found at iter:     $(stats.iter_found_solution)")
    end

    if stats.infeasible_exit
        println("Exit reason:       infeasible FW solution")
    end
    
    println("="^80 * "\n")

    return (SCIP.SCIP_OKAY, result)
end

function submit_solution_to_scip(
    scip::Ptr{SCIP.SCIP_},
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
    lp_cols::Vector{Ptr{SCIP.SCIP_COL}},
    solution::Vector{Float64},
    nvars::Int32
)::Bool

    sol_ptr = Ref{Ptr{SCIP.SCIP_SOL}}()
    SCIP.SCIPcreateSol(scip, sol_ptr, heur_ptr)
    sol = sol_ptr[]
    
    # Set solution values
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)
        SCIP.SCIPsetSolVal(scip, sol, var, solution[j])
    end
    
    # Try to add solution
    stored = Ref{SCIP.SCIP_Bool}()
    SCIP.SCIPtrySol(scip, sol, SCIP.TRUE, SCIP.FALSE, SCIP.TRUE, SCIP.TRUE, SCIP.TRUE, stored)
    
    if stored[] == SCIP.TRUE
        println("Solution accepted!")
        return true
    else
        SCIP.SCIPfreeSol(scip, sol_ptr)
        return false
    end
end