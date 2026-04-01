# Helper function to assign Frank-Wolfe variant.
# All variants share the same return signature (x, v, primal, dual_gap, traj).
# Note: :away, :blended_pairwise, and :blended maintain an active set internally and converge faster on smooth objectives; 
# :vanilla is the safest choice for non-smooth (manhattan) objectives.
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

    # DEBUG: Print initial LP solution details
    if DEBUG_VERBOSE
        println("[DEBUG] Initial LP solution:")
        for j in 1:nvars
            var = SCIP.SCIPcolGetVar(lp_cols[j])
            @printf("  x[%d] = %.3f\n", j, current_solution[j])
        end
    end

    lp_obj = 0.0
    for j in 1:nvars
        lp_obj += current_solution[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j]))
    end

    println("Initial LP objective: $(round(lp_obj, digits=4))")

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

    else
        error("Unknown projection norm: $(heur.projection_norm)")
    end

    # Hash function for cycle detection (only hashes binary variable values)
    function hash_rounded(x_round, binary_indices)
        return hash(tuple((x_round[i] for i in binary_indices)...))
    end

    # Rounding function using custom threshold
    function round_with_threshold(x, binary_indices, threshold)
        x_round = copy(x)
        for i in binary_indices
            x_round[i] = x[i] >= threshold ? 1.0 : 0.0
        end
        return x_round
    end

    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(current_solution)  # Initial LP feasible solution
    prev_x = copy(x)  # To calculate distance moved
    last_feasible_x = copy(current_solution)  # Always keep a feasible point for restarts

    # Cycle detection: store hashes of all visited rounded solutions
    rounded_cache = Set{UInt64}()
    restart_count = 0

    # Store first feasible solution found
    found_solution = nothing

    for iter in 1:DEF_FP_MAX_ITER
        println("\n--- FPFW Iteration $iter ---")
        stats.fp_iterations = iter

        # Check time limit
        elapsed = time() - total_start_time
        if elapsed > DEF_FW_TIME_LIMIT
            println("FW time limit reached ($(round(elapsed, digits=2))s > $(DEF_FW_TIME_LIMIT)s)")
            break
        end

        # Step 1: Rounding LP feasible solution with custom threshold
        x_round = round_with_threshold(x, binary, heur.rounding_threshold)

        # DEBUG: Print before/after rounding
        if DEBUG_VERBOSE
            println("[DEBUG] Rounding (threshold=$(heur.rounding_threshold)):")
            for i in binary
                @printf("  x[%d]: %.3f -> %.1f\n", i, x[i], x_round[i])
            end
        end

        # Cycle detection: check if we've visited this rounded solution before
        h = hash_rounded(x_round, binary)
        if h in rounded_cache
            restart_count += 1
            println("Cycle detected at iteration $iter (restart #$restart_count)")

            if restart_count >= DEF_MAX_RESTARTS
                println("Maximum restarts ($DEF_MAX_RESTARTS) reached, stopping.")
                break
            end

            # Perturb the rounding target directly (flip some binary values)
            for i in binary
                if rand() < DEF_PERTURB_FRACTION
                    x_round[i] = 1.0 - x_round[i]  # Flip 0-1
                end
            end
            # Reset x to last known feasible point
            x .= last_feasible_x
            println("Perturbed rounding target, restarting from feasible point...")
            # Don't skip - proceed with FW projection toward the perturbed target
        end
        push!(rounded_cache, h)

        fw_traj = Vector{Any}()
        fw_callback = FrankWolfe.make_trajectory_callback(nothing, fw_traj)

        # Step 2: Projection using Frank-Wolfe
        ls = heur.line_search == :backtracking ? FrankWolfe.Backtracking() : FrankWolfe.Agnostic()
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
            epsilon = 1e-7,
            callback = fw_callback
        )
        stats.fw_time += time() - fw_start

        fw_iters = length(fw_traj)
        stats.fw_iterations += fw_iters

        # DEBUG: Print FW projection result
        if DEBUG_VERBOSE
            total_dist = sum(abs(x_new[i] - x_round[i]) for i in binary)
            @printf("[DEBUG] FW done: %d iters, dist to rounding target = %.6f\n", fw_iters, total_dist)
        end

        # Step 3: Check feasibility and integrality and calculate current objective value
        is_integral = check_integrality(x_new, binary)
        is_feasible = check_feasibility(scip, x_new, col_to_idx)
        distance = sqrt(sum((x_new[i] - prev_x[i])^2 for i in binary))

        obj_val = 0.0
        for j in 1:nvars
            obj_val += x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j]))
        end

        # DEBUG: Print integrality check details
        if DEBUG_VERBOSE
            @printf("[DEBUG] integral=%s, feasible=%s, obj=%.6f, dist=%.6f\n",
                is_integral, is_feasible, obj_val, distance)
        end

        # Just for debugging
        if !is_feasible
            println("FrankWolfe return infeasible solution!")
            break
        end

        if is_integral
            stats.solution_found = true
            stats.final_objective = obj_val

            found_solution = copy(x_new)

            # DEBUG: Print the found solution in detail
            if DEBUG_VERBOSE
                println("[DEBUG] Solution values:")
                for j in 1:nvars
                    var = SCIP.SCIPcolGetVar(lp_cols[j])
                    @printf("  x[%d] = %.6f (obj=%.4f)\n", j, x_new[j], SCIP.SCIPvarGetObj(var))
                end
            end

            # Submit solution to SCIP
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_new, nvars)
                result = SCIP.SCIP_FOUNDSOL
                break  # Exit after finding a feasible solution
            else
                println("Solution rejected by SCIP")
                break
            end
        else
            # Fixed-point detection: solution barely changed
            if distance < DEF_TOLERANCE
                restart_count += 1
                println("Converged to fixed point at iteration $iter (restart #$restart_count)")

                if restart_count >= DEF_MAX_RESTARTS
                    println("Maximum restarts ($DEF_MAX_RESTARTS) reached, stopping.")
                    break
                end

                # Reset to feasible point and let next iteration create a new rounding target
                x .= last_feasible_x
                # Slightly perturb x within feasible bounds to get different rounding
                for i in binary
                    x[i] = clamp(x[i] + 0.1 * (rand() - 0.5), 0.0, 1.0)
                end
                println("Perturbed and reset, continuing...")
                continue
            end

            # Continue with the projected solution for next FW iteration
            prev_x .= x
            x .= x_new
            last_feasible_x .= x_new  # Store as last known feasible point
        end
    end

    # Print final summary
    total_time = time() - total_start_time
    stats.total_time = total_time
    println("\n" * "="^80)
    println("FPFW HEURISTIC SUMMARY")
    println("="^80)
    println("Total time:        $(round(total_time, digits=2)) seconds")
    println("FP iterations:     $(stats.fp_iterations)")
    println("FW iterations:     $(stats.fw_iterations)")
    println("FW time:           $(round(stats.fw_time, digits=2)) seconds")
    println("Restarts (cycles): $restart_count")
    println("Solution found:    $(stats.solution_found)")

    if found_solution !== nothing
        println("\nFIRST SOLUTION FOUND:")
        println("  Objective value: $(stats.final_objective)")
        println("  Solution vector (first 20 vars): $(found_solution[1:min(20, length(found_solution))])")
    else
        println("\nNo feasible integer solution found.")
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