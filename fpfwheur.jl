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
    stats = FPFWStats()
    startTime = time()

    # SCIP is initially given DEF_SCIP_TIME_LIMIT (300s) to solve the LP relaxation.                                                                                                                                  
    # Once the heuristic is called, extend to DEF_GLOBAL_TIME_LIMIT (480s) so SCIP doesn't terminate while the heuristic is running.
    SCIP.SCIPsetRealParam(scip, "limits/time", DEF_GLOBAL_TIME_LIMIT)

    # Extract LP data
    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    lp_cols, col_to_idx, binary, integer, current_solution = extract_lp_data(scip, nvars)
    intIndices = [binary; integer]

    # Build LMO from current LP
    # Currently not useful, since we set the heuristic to only run once
    if heur.lmo === nothing
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
    end

    # Log initial LP solve info
    stats.dualBound = SCIP.SCIPgetDualbound(scip)
    lpRootIter = SCIP.SCIPgetNRootLPIterations(scip)
    nFracVars = count(i -> abs(current_solution[i] - round(current_solution[i])) > DEF_TOLERANCE, intIndices)
    rootTime = SCIP.SCIPgetSolvingTime(scip)

    printstyled("[initialSolve]\n", color=:cyan)
    @printf("Initial LP: lpiter=%d obj=%.2f frac=%d/%d time=%.2fs\n", 
            lpRootIter, stats.dualBound, nFracVars, length(intIndices), rootTime)

    if DEBUG_VERBOSE
        printstyled("[debug info]\n", color=:yellow)
        for j in 1:nvars
            @printf("  x[%d] = %.3f\n", j, current_solution[j])
        end        
    end

    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(current_solution)       # Initial LP feasible solution
    prev_x = copy(current_solution)  # To calculate distance moved
    x_round = zeros(Float64, nvars)  # Preallocated buffer for rounding target
    x_after = zeros(Float64, nvars)  # Preallocated buffer for post-step x in callback
    x_temp  = zeros(Float64, nvars)  # Preallocated buffer for randomized rounding

    # Active set for warm-starting away/blended variants
    activeSet = nothing             # Preallocate active set for warm-starting

    # Stagnation detection
    bestIntGap = Inf
    stagnationCount = 0

    # Randomized rounding parameters
    attempts = max(1, length(intIndices))

    # FW escape flag and buffer
    fwEscaped = false
    x_round_escape = zeros(Float64, nvars)

    # TODO: store the best solution found across iterations, not just the first one
    found_solution = nothing

    f, grad! = build_fw_functions(heur.config.projectionNorm, intIndices)
    ls = build_line_search(heur.config.lineSearch)

    # FW step trajectory within one FP iteration: (step, x, objective)
    fwTraj = Vector{Tuple{Int, Vector{Float64}, Float64}}()
    fwCallback = (state, args...) -> begin
        # Skip FW bookkeeping steps where d and gamma are not meaningful
        if state.step_type === FrankWolfe.ST_LAST || state.step_type === FrankWolfe.ST_POSTPROCESS
            return true
        end

        if state.d === nothing || state.gamma === nothing
            return true
        end

        # Compute next iterate and log it
        x_after .= state.x .- state.gamma .* state.d
        push!(fwTraj, (state.t, copy(x_after), state.primal))

        # Check if x_after rounds to a different target than x_round
        # If so, stop FW early and use x_after as the new starting point
        if DEF_FW_ESCAPE
            round_solution!(x_round_escape, x_after, intIndices, DEF_ROUNDING_THRESHOLD)
            if !are_equal_vectors(intIndices, x_round_escape, x_round)
                fwEscaped = true
                return false  # stop FW iter early
            end
        end

        return true  # continue FW iter
    end

    # Main FPFW loop
    while true
        stats.fpIterations += 1
        if DEBUG_VERBOSE
            printstyled("\n--- FPFW Iteration $(stats.fpIterations) ---\n"; color=:blue)
        end

        # Check time limit
        if time() - heur.globalStartTime > DEF_GLOBAL_TIME_LIMIT
            stats.exitReason = :time_limit
            break
        end

        if heur.config.randRound && stats.fpIterations > 1  # skip randomized rounding in the first iteration to save time
            rrStartTime = time()
            for _ in 1:attempts
                if time() - rrStartTime > DEF_RR_TIME_LIMIT
                    break
                end

                for i in intIndices
                    frac = x[i] - floor(x[i])
                    x_temp[i] = rand() < frac ? ceil(x[i]) : floor(x[i])
                end

                if check_feasibility(scip, x_temp, col_to_idx)
                    if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_temp, nvars)
                        stats.solutionFound = true
                        stats.exitReason = :rr_solution_found
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

            stats.rrTime += time() - rrStartTime
                
            if stats.solutionFound
                break
            end
        end

        if stagnationCount >= DEF_MAX_STAGNATION
            stats.restartCount += 1

            if stats.restartCount >= DEF_MAX_RESTARTS
                stats.exitReason = :restart_limit
                break
            end

            if DEBUG_VERBOSE
                println("Cycle detected at iteration $(stats.fpIterations) (restart #$(stats.restartCount))")
                println("Perturbing:")
            end

            Random.seed!(DEF_RANDOM_SEED + stats.restartCount)  # change seed each restart for reproducibility
            perturb_solution!(x, x_round, binary, integer, lp_cols)

            stagnationCount = 0
            bestIntGap = Inf
        else
            round_solution!(x_round, x, intIndices, DEF_ROUNDING_THRESHOLD)
        end

        # Check if rounded solution is feasible
        if check_feasibility(scip, x_round, col_to_idx)
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_round, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
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
            for i in intIndices
                @printf("  x[%d]: %.3f -> %.1f\n", i, x[i], x_round[i])
            end
            println("------------------------\n")
        end

        # Reset FW escape flag and trajectory
        fwEscaped = false
        empty!(fwTraj)

        # Step 2: "Projection" using Frank-Wolfe
        remainingTime = DEF_GLOBAL_TIME_LIMIT - (time() - heur.globalStartTime)
        fwStartTime = time()

        fwResult = run_fw(
            heur.config.fwVariant,
            x -> f(x, x_round),
            (storage, x) -> grad!(storage, x, x_round),
            heur.lmo,
            x,
            activeSet,
            heur.config.warmStart,
            ls,
            fwCallback,
            remainingTime
        )

        stats.fwTime += time() - fwStartTime
        fwIters = length(fwTraj)
        stats.fwIterations += fwIters

        # If FW gives us intermediate solutions via callback that escape the rounding target, we immediately jump to next iteration with the escaped solution
        if fwEscaped
            x .= x_after
            prev_x .= x_after
            continue  # skip rest of checks and go to next FPFW iteration
        else
            x_new = fwResult.x
            if heur.config.warmStart && heur.config.fwVariant !== :vanilla
                activeSet = fwResult.activeSet
            end
        end

        # Cycle detection: check if we've visited this solution before
        intGap = f(x_new, x_round)

        if is_lower_than(intGap, bestIntGap)
            bestIntGap = intGap
            stagnationCount = 0  # reset stagnation count if we made progress
        else
            stagnationCount += 1
        end

        if DEBUG_VERBOSE
            for (t, xk, fk) in fwTraj
                println("--- FW Step $t ---")
                println("Solution:")
                for i in intIndices
                    @printf("  x[%d]: %.3f\n", i, xk[i])
                end
                @printf("Objective = %.3f\n\n", fk)
            end
        end

        # Step 3: Check feasibility, integrality, distance moved, and objective value
        isIntegral = check_integrality(x_new, intIndices)
        isFeasible = check_feasibility(scip, x_new, col_to_idx)

        if DEBUG_VERBOSE
            obj_val = sum(x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
            dist_moved = f(x_new, prev_x)
            int_gap = f(x_new, x_round)
            nfrac_binary = count(i -> abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE, binary)
            nfrac_integer = count(i -> abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE, integer)

            @printf("[Iter %d] obj=%.3f | distance_moved=%.3f | integrality_gap=%.3f | frac_bin=%d/%d | frac_int=%d/%d | fw_iters=%d | integral=%s |feasible=%s\n",
                stats.fpIterations, obj_val, dist_moved, int_gap, nfrac_binary, length(binary), nfrac_integer, length(integer), fwIters, isIntegral, isFeasible)
        end

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isFeasible
            stats.exitReason = :infeasible_fw
            break
        end

        if isIntegral
            # Submit solution to SCIP
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_new, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
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
                stats.exitReason = :solution_rejected
                break
            end
        else
            # Continue with the projected solution for next FW iteration
            prev_x .= x_new
            x .= x_new
        end
    end

    # Print final summary
    stats.heurTime = time() - startTime
    stats.totalTime = time() - heur.globalStartTime
    stats.primalBound = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    print_heuristic_summary(stats)

    return (SCIP.SCIP_OKAY, result)
end
