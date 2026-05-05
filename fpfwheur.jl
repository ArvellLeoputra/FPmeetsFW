# Main FPFW Heuristic Implementation
function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}

    @assert SCIP.SCIPhasCurrentNodeLP(scip) == SCIP.TRUE  # always true, since we set timing to DURINGLPLOOP
    
    heur.called += 1
    if heur.called > 1
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end
    
    result = SCIP.SCIP_DIDNOTFIND

    # SCIP is initially given DEF_SCIP_TIME_LIMIT (300s) to solve the LP relaxation.                                                                                                                                  
    # Once the heuristic is called, extend to DEF_GLOBAL_TIME_LIMIT (480s) so SCIP doesn't terminate while the heuristic is running.
    SCIP.SCIPsetRealParam(scip, "limits/time", DEF_GLOBAL_TIME_LIMIT)

    stats = FPFWStats()
    heur_start_time = time()

    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)

    # Build LMO from current LP
    # Currently not useful, since we set the heuristic to only run once
    if heur.lmo === nothing
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
    end

    lp_cols, col_to_idx, binary, integer, current_solution = extract_lp_data(scip, nvars)
    all_integers = [binary; integer]

    if DEBUG_VERBOSE
        println("Initial solution:")
        for j in 1:nvars
            @printf("  x[%d] = %.3f\n", j, current_solution[j])
        end
    end

    lp_objective = sum(current_solution[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j]))
                             for j in 1:nvars)

    if DEBUG_VERBOSE
        println("Initial LP objective: $(round(lp_objective, digits=4))")

        nfrac_binary = count(i -> abs(current_solution[i] - round(current_solution[i])) > DEF_TOLERANCE, binary)
        println("Non-integral binary vars: $nfrac_binary/$(length(binary))")

        nfrac_integer = count(i -> abs(current_solution[i] - round(current_solution[i])) > DEF_TOLERANCE, integer)
        println("Non-integral integer vars: $nfrac_integer/$(length(integer))")

    end

    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(current_solution)       # Initial LP feasible solution
    prev_x = copy(current_solution)  # To calculate distance moved
    x_round = zeros(Float64, nvars)  # Preallocated buffer for rounding target
    x_after = zeros(Float64, nvars)  # Preallocated buffer for post-step x in callback
    x_temp  = zeros(Float64, nvars)  # Preallocated buffer for randomized rounding

    # Active set for warm-starting away/blended variants
    active_set = nothing             # Preallocate active set for warm-starting

    # Stagnation detection
    best_int_gap = Inf
    stagnation_count = 0

    # FW escape flag and buffer
    fw_escaped = false
    x_round_escape = zeros(Float64, nvars)

    # TODO: store the best solution found across iterations, not just the first one
    found_solution = nothing

    f, grad! = build_fw_functions(heur.config.projection_norm, all_integers)
    ls = build_line_search(heur.config.line_search)

    # FW step trajectory within one FP iteration: (step, x, objective)
    fw_traj = Vector{Tuple{Int, Vector{Float64}, Float64}}()
    fw_callback = (state, args...) -> begin
        # Skip FW bookkeeping steps where d and gamma are not meaningful
        if state.step_type === FrankWolfe.ST_LAST || state.step_type === FrankWolfe.ST_POSTPROCESS
            return true
        end

        if state.d === nothing || state.gamma === nothing
            return true
        end

        # Compute next iterate and log it
        x_after .= state.x .- state.gamma .* state.d
        push!(fw_traj, (state.t, copy(x_after), state.primal))

        # Check if x_after rounds to a different target than x_round
        # If so, stop FW early and use x_after as the new starting point
        round_solution!(x_round_escape, x_after, all_integers, DEF_ROUNDING_THRESHOLD)
        if !are_equal_vectors(all_integers, x_round_escape, x_round)
            fw_escaped = true
            return false  # stop FW iter early
        end

        return true  # continue FW iter
    end

    attempts = max(1, length(all_integers))

    iter = 0
    while true
        iter += 1
        if DEBUG_VERBOSE
            println("\n--- FPFW Iteration $iter ---")
        end

        # Check time limit
        if time() - heur.global_start_time > DEF_GLOBAL_TIME_LIMIT
            stats.exit_reason = :time_limit
            break
        end

        if heur.config.rand_round && iter > 1  # skip randomized rounding in the first iteration to save time
            rr_iter_start = time()
            for _ in 1:attempts
                if time() - rr_iter_start > DEF_RR_TIME_LIMIT
                    break
                end

                for i in all_integers
                    frac = x[i] - floor(x[i])
                    x_temp[i] = rand() < frac ? ceil(x[i]) : floor(x[i])
                end

                if check_feasibility(scip, x_temp, col_to_idx)
                    if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_temp, nvars)
                        stats.solution_found = true
                        stats.exit_reason = :rr_solution_found
                        stats.iter_found_solution = iter
                        result = SCIP.SCIP_FOUNDSOL

                        if DEBUG_VERBOSE
                            println()
                            println("End solution:")
                            for j in 1:nvars
                                @printf("  x[%d] = %.3f\n", j, x_temp[j])
                            end
                        end
                        break
                    end
                end
            end

            stats.rr_time += time() - rr_iter_start
                
            if stats.solution_found
                break
            end
        end

        stats.fp_iterations = iter

        # Step 1: Rounding LP feasible solution with custom threshold
        round_solution!(x_round, x, all_integers, DEF_ROUNDING_THRESHOLD)

        # Check if rounded solution is feasible
        if check_feasibility(scip, x_round, col_to_idx)
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_round, nvars)
                stats.solution_found = true
                stats.exit_reason = :solution_found
                stats.iter_found_solution = iter
                result = SCIP.SCIP_FOUNDSOL

                if DEBUG_VERBOSE
                    println()
                    println("End solution:")
                    for j in 1:nvars
                        @printf("  x[%d] = %.3f\n", j, x_round[j])
                    end
                end
                break
            end
        end

        if DEBUG_VERBOSE
            println("Rounding(threshold=$(DEF_ROUNDING_THRESHOLD)):")
            for i in all_integers
                @printf("  x[%d]: %.3f -> %.1f\n", i, x[i], x_round[i])
            end
            println("------------------------\n")
        end

        fw_escaped = false
        empty!(fw_traj)

        # Step 2: Projection using Frank-Wolfe
        remaining_time = DEF_GLOBAL_TIME_LIMIT - (time() - heur.global_start_time)
        fw_start = time()

        fw_result = run_fw(
            heur.config.fw_variant,
            x -> f(x, x_round),
            (storage, x) -> grad!(storage, x, x_round),
            heur.lmo,
            x,
            active_set,
            heur.config.warm_start,
            ls,
            fw_callback,
            remaining_time
        )

        # If FW gives us intermediate solutions via callback that escape the rounding target, we immediately jump to next iteration with the escaped solution
        if fw_escaped
            x .= x_after
            prev_x .= x_after
            continue
        else
            x_new = fw_result.x
            if heur.config.warm_start && heur.config.fw_variant !== :vanilla
                active_set = fw_result.active_set
            end
        end

        # Cycle detection: check if we've visited this solution before
        int_gap = f(x_new, x_round)
        if is_lower_than(int_gap, best_int_gap)
            best_int_gap = int_gap
            stagnation_count = 0
        else
            stagnation_count += 1
        end

        if stagnation_count >= DEF_MAX_STAGNATION
            stats.restarts += 1

            if stats.restarts >= DEF_MAX_RESTARTS
                stats.exit_reason = :restart_limit
                break
            end

            if DEBUG_VERBOSE
                println("Cycle detected at iteration $iter (restart #$(stats.restarts))")
                println("Perturbing:")
            end

            for i in binary
                if rand() < DEF_PERTURB_FRACTION
                    x_round[i] = 1.0 - x_round[i]
                end
            end

            for i in integer
                if rand() < DEF_PERTURB_FRACTION
                    var = SCIP.SCIPcolGetVar(lp_cols[i])
                    lb = SCIP.SCIPvarGetLbLocal(var)
                    ub = SCIP.SCIPvarGetUbLocal(var)
                    r = rand()

                    newval = if (ub - lb) < DEF_BIGBIGM
                        floor(lb + (1 + ub - lb) * r)
                    elseif (x_round[i] - lb) < DEF_BIGM
                        lb + (2 * DEF_BIGM - 1) * r
                    elseif (ub - x_round[i]) < DEF_BIGM
                        ub - (2 * DEF_BIGM - 1) * r
                    else
                        x[i] + (2 * DEF_BIGM - 1) * r - DEF_BIGM
                    end

                    x_round[i] = clamp(floor(newval), lb, ub)
                end
            end

            if DEBUG_VERBOSE
                for i in binary
                    @printf("  bin x[%d]: %.1f\n", i, x_round[i])
                end

                for i in integer
                    @printf("  int x[%d]: %.1f\n", i, x_round[i])
                end
                println()
            end

            stagnation_count = 0
            best_int_gap = Inf

            # Re-run FW with perturbed x_round so perturbation takes effect immediately
            remaining_time = DEF_GLOBAL_TIME_LIMIT - (time() - heur.global_start_time)
            fw_escaped = false
            empty!(fw_traj)

            fw_result_perturbed = run_fw(
                heur.config.fw_variant,
                x -> f(x, x_round),
                (storage, x) -> grad!(storage, x, x_round),
                heur.lmo,
                x,
                active_set,
                heur.config.warm_start,
                ls,
                fw_callback,
                remaining_time
            )

            if fw_escaped
                x .= x_after
            else
                x .= fw_result_perturbed.x
                if heur.config.warm_start && heur.config.fw_variant !== :vanilla
                    active_set = fw_result_perturbed.active_set
                end
            end

            prev_x .= x
        end

        stats.fw_time += time() - fw_start
        fw_iters = length(fw_traj)
        stats.fw_iterations += fw_iters

        # Check time limit after FW returns
        if time() - heur.global_start_time > DEF_GLOBAL_TIME_LIMIT
            stats.exit_reason = :time_limit
            break
        end

        if DEBUG_VERBOSE
            for (t, xk, fk) in fw_traj
                println("--- FW Step $t ---")
                println("Solution:")
                for i in all_integers
                    @printf("  x[%d]: %.3f\n", i, xk[i])
                end
                @printf("Objective = %.3f\n\n", fk)
            end
        end

        # Step 3: Check feasibility, integrality, distance moved, and objective value
        is_integral = check_integrality(x_new, all_integers)
        is_feasible = check_feasibility(scip, x_new, col_to_idx)
        obj_val = sum(x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)

        if DEBUG_VERBOSE
            nfrac_binary = count(i -> abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE, binary)
            nfrac_integer = count(i -> abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE, integer)
            dist_moved = f(x_new, prev_x)
            int_gap = f(x_new, x_round)

            @printf("[Iter %d] obj=%.3f | distance_moved=%.3f | integrality_gap=%.3f | frac_bin=%d/%d | frac_int=%d/%d | fw_iters=%d | integral=%s |feasible=%s\n",
                iter, obj_val, dist_moved, int_gap, nfrac_binary, length(binary), nfrac_integer, length(integer), fw_iters, is_integral, is_feasible)
        end

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !is_feasible
            stats.exit_reason = :infeasible_fw
            break
        end

        if is_integral
            # Submit solution to SCIP
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_new, nvars)
                stats.solution_found = true
                stats.exit_reason = :solution_found
                stats.iter_found_solution = iter
                result = SCIP.SCIP_FOUNDSOL

                if DEBUG_VERBOSE
                    println()
                    println("End solution:")
                    for j in 1:nvars
                        @printf("  x[%d] = %.3f\n", j, x_new[j])
                    end
                end
                break
            else
                stats.exit_reason = :solution_rejected
                break
            end
        else
            # Continue with the projected solution for next FW iteration
            prev_x .= x_new
            x .= x_new
        end
    end

    # Print final summary
    stats.heur_time = time() - heur_start_time
    total_time = time() - heur.global_start_time
    objective = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
    gap = Float64(SCIP.SCIPgetGap(scip))

    print_heuristic_summary(stats, total_time, objective, gap)

    return (SCIP.SCIP_OKAY, result)
end
