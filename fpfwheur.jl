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
    lp_cols, lp_rows, col_to_idx, binary, integer, current_solution = extract_lp_data(scip, nvars, nrows)
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
        println("Initial LP solution:")
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
    restarted = false

    # Randomized rounding parameters
    attempts = min(DEF_RAND_FEAS_ITER_LIMIT, length(intIndices))

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
            round_solution!(x_round_escape, x_after, intIndices, heur.config.randRound)
            if !are_equal_vectors(intIndices, x_round_escape, x_round)
                fwEscaped = true    
                
                obj = sum(x_after[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
                step = f(x_after, prev_x)
                intGap = f(x_after, x_round)
                nFrac = count(i -> abs(x_after[i] - round(x_after[i])) > DEF_TOLERANCE, intIndices)
                escFwIters = state.t
                iterTime = time() - iterStartTime

                if DEBUG_VERBOSE
                    @printf("FW escaped rounding target at step %d: obj=%.4f projObj=%.4f step=%.4f nFrac=%d\n",
                            escFwIters, obj, intGap, step, nFrac)
                else
                    print_row!(pump_display, stats.pumpIterations, obj, intGap, step, nFrac, escFwIters, iterTime, "escape")
                end

                return false  # stop FW iter early
            end
        end

        if DEF_FW_MAX_ITER > 0 && state.t >= DEF_FW_MAX_ITER                                                                                                                                                                               
            return false                                                                                                                                                                                                                   
        end
        
        return true  # continue FW iter
    end

    if !DEBUG_VERBOSE
        printstyled("[pump]\n", color=:cyan)
    end

    pump_display = PumpDisplay(PumpDisplayColumn[])                                                                                                                                                                                   
    add_column!(pump_display, "pumpIter", 10)                                                                                                                                                                                         
    add_column!(pump_display, "obj", 15, 4)                                                                                                                                                                                           
    add_column!(pump_display, "projObj", 15, 4)                                                                                                                                                                                       
    add_column!(pump_display, "step", 15, 4)                                                                                                                                                                                          
    add_column!(pump_display, "nFrac", 8)                                                                                                                                                                                             
    add_column!(pump_display, "fwIters", 10)                                                                                                                                                                                          
    add_column!(pump_display, "time", 10, 2)                                                                                                                                                                                          
    add_column!(pump_display, "flag", 18)   
                                                                                                                                                                                        
    if !DEBUG_VERBOSE
        print_header!(pump_display)
    end

    # Main FPFW loop
    while true
        stats.pumpIterations += 1
        
        if DEBUG_VERBOSE
            printstyled("\nFPFW Iteration $(stats.pumpIterations)\n"; color=:blue)
        end

        iterStartTime = time()
        restarted = false

        # Check time limit
        if iterStartTime - heur.globalStartTime > DEF_GLOBAL_TIME_LIMIT
            stats.exitReason = :time_limit
            break
        end

        if heur.config.randFeasCheck && stats.pumpIterations > 1  # skip randomized rounding in the first iteration to save time
            rrStartTime = time()
            for _ in 1:attempts
                for i in intIndices
                    frac = x[i] - floor(x[i])
                    x_temp[i] = rand() < frac ? ceil(x[i]) : floor(x[i])
                end

                if check_feasibility(scip, lp_rows, lp_cols, x_temp, col_to_idx)
                    if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_temp, nvars)
                        stats.solutionFound = true
                        stats.exitReason = :rr_solution_found
                        result = SCIP.SCIP_FOUNDSOL

                        if DEBUG_VERBOSE
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
                obj = sum(x_temp[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)                                                                                                                        
                iterTime = time() - iterStartTime                                                                                                                                                                                 
                if !DEBUG_VERBOSE
                    print_row!(pump_display, stats.pumpIterations, obj, NaN, NaN, 0, 0, iterTime, "randFeasCheck")
                end
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
                println("Cycle detected at iteration $(stats.pumpIterations) (restart #$(stats.restartCount))")
                println("Perturbing:")
            end

            Random.seed!(DEF_RANDOM_SEED + stats.restartCount)  # change seed each restart for reproducibility
            perturb_solution!(x, x_round, binary, integer, lp_cols)

            stagnationCount = 0
            bestIntGap = Inf
            restarted = true
        else
            round_solution!(x_round, x, intIndices, heur.config.randRound)
        end

        if DEBUG_VERBOSE
            println("Rounding:")
            for i in intIndices
                @printf("  x[%d]: %.3f -> %d\n", i, x[i], Int(x_round[i]))
            end
        end

        # Check if rounded solution is feasible
        if check_feasibility(scip, lp_rows, lp_cols, x_round, col_to_idx)
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_round, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
                result = SCIP.SCIP_FOUNDSOL
                obj = sum(x_round[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
                iterTime = time() - iterStartTime
                if !DEBUG_VERBOSE
                    print_row!(pump_display, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restart+feasRound" : "feasRound")
                end
                break
            end
        end

        feasible, sol = lp_diving!(scip, lp_cols, intIndices, x_round, nvars)                                                                                                                                                              
        if feasible
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, sol, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
                result = SCIP.SCIP_FOUNDSOL
                obj = sum(sol[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
                iterTime = time() - iterStartTime
                if !DEBUG_VERBOSE
                    print_row!(pump_display, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restart+divingLP" : "divingLP")
                end
                break
            end
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
            x_new = isempty(fwTraj) ? fwResult.x : fwTraj[end][2]
            if heur.config.warmStart && heur.config.fwVariant !== :vanilla
                activeSet = fwResult.active_set
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
                t > DEF_FW_MAX_ITER && DEF_FW_MAX_ITER > 0 && continue
                println("--- FW Step $t ---")
                println("Solution:")
                for i in intIndices
                    @printf("  x[%d]: %.3f\n", i, xk[i])
                end
                @printf("Objective = %.3f\n", fk)
            end
        end

        # Step 3: Check feasibility, integrality, distance moved, and objective value
        isIntegral = check_integrality(x_new, intIndices)
        isFeasible = check_feasibility(scip, lp_rows, lp_cols, x_new, col_to_idx)

        obj = sum(x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j])) for j in 1:nvars)
        step = f(x_new, prev_x)
        nFrac = count(i -> abs(x_new[i] - round(x_new[i])) > DEF_TOLERANCE, intIndices)
        iterTime = time() - iterStartTime

        if !DEBUG_VERBOSE
            print_row!(pump_display, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, restarted ? "restart" : "")
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
