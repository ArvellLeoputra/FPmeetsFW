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

    # Get LP data
    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    lpCols, lpRows, colDict, binIdx, gIntIdx, initSol = getLPData(scip, nvars, nrows)
    intIdx = [binIdx; gIntIdx]

    # Build LMO from current LP
    # Currently not useful, since we set the heuristic to only run once
    if heur.lmo === nothing
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
    end

    # Log initial LP solve info
    stats.dualBound = SCIP.SCIPgetDualbound(scip)
    lpRootIter = SCIP.SCIPgetNRootLPIterations(scip)
    nFracVars = countFracVars(intIdx, initSol)
    rootTime = SCIP.SCIPgetSolvingTime(scip)

    printstyled("[initialSolve]\n", color=:cyan)
    @printf("Initial LP: lpiter=%d obj=%.2f frac=%d/%d time=%.2fs\n", 
            lpRootIter, stats.dualBound, nFracVars, length(intIdx), rootTime)

    if DEBUG_VERBOSE
        printstyled("[debug info]\n", color=:yellow)
        println("Initial LP solution:")
        for j in 1:nvars
            @printf("  x[%d] = %.3f\n", j, initSol[j])
        end        
    end

    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(initSol)       # Initial LP feasible solution
    xPrev = copy(initSol)   # To calculate distance moved
    xRound = zeros(Float64, nvars)  # Preallocated buffer for rounding target
    xAfter = zeros(Float64, nvars)  # Preallocated buffer for post-step x in callback
    xTemp  = zeros(Float64, nvars)  # Preallocated buffer for randomized rounding

    # Active set for warm-starting away/blended variants
    activeSet = nothing             # Preallocate active set for warm-starting

    # Stagnation detection
    bestIntGap = Inf
    stagnationCount = 0
    restarted = false
    visitedRounded = Set{UInt}()

    # Randomized rounding parameters
    attempts = min(DEF_RAND_FEAS_ITER_LIMIT, length(intIdx))

    # FW escape flag and buffer
    fwEscaped = false
    xRoundEscape = zeros(Float64, nvars)

    # TODO: store the best solution found across iterations, not just the first one
    foundSolution = nothing

    f, grad! = buildFWFunctions(heur.config.projNorm, intIdx)
    ls = buildLineSearch(heur.config.lineSearch)

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
        xAfter .= state.x .- state.gamma .* state.d
        push!(fwTraj, (state.t, copy(xAfter), state.primal))

        # Check if xAfter rounds to a different target than xRound
        # If so, stop FW early and use xAfter as the new starting point
        if DEF_FW_ESCAPE
            roundSolution!(xRoundEscape, xAfter, intIdx, heur.config.randRound)
            if areSolutionsEqual(intIdx, xRoundEscape, xRound)
                fwEscaped = true    
                
                obj = sum(xAfter[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:nvars)
                step = f(xAfter, xPrev)
                intGap = f(xAfter, xRound)
                nFrac = count(i -> abs(xAfter[i] - round(xAfter[i])) > DEF_INT_TOLERANCE, intIdx)
                escFwIters = state.t
                iterTime = time() - iterStartTime

                if DEBUG_VERBOSE
                    @printf("FW escaped rounding target at step %d: obj=%.4f projObj=%.4f step=%.4f nFrac=%d\n",
                            escFwIters, obj, intGap, step, nFrac)
                else
                    printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, escFwIters, iterTime, "escape")
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

    pumpDisplay = PumpDisplay(PumpDisplayColumn[])                                                                                                                                                                                   
    addColumn!(pumpDisplay, "pumpIter", 10)                                                                                                                                                                                         
    addColumn!(pumpDisplay, "obj", 15, 4)                                                                                                                                                                                           
    addColumn!(pumpDisplay, "projObj", 15, 4)                                                                                                                                                                                       
    addColumn!(pumpDisplay, "step", 15, 4)                                                                                                                                                                                          
    addColumn!(pumpDisplay, "nFrac", 8)                                                                                                                                                                                             
    addColumn!(pumpDisplay, "fwIters", 10)                                                                                                                                                                                          
    addColumn!(pumpDisplay, "time", 10, 2)                                                                                                                                                                                          
    addColumn!(pumpDisplay, "flag", 18)   
                                                                                                                                                                                        
    if !DEBUG_VERBOSE
        printHeader!(pumpDisplay)
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
                for i in intIdx
                    frac = x[i] - floor(x[i])
                    xTemp[i] = rand() < frac ? ceil(x[i]) : floor(x[i])
                end

                if submitSolution(scip, heur_ptr, lpCols, xTemp, nvars)
                    stats.solutionFound = true
                    stats.exitReason = :rr_solution_found
                    result = SCIP.SCIP_FOUNDSOL

                    if DEBUG_VERBOSE
                        println("End solution:")
                        for j in 1:nvars
                            @printf("  x[%d] = %.3f\n", j, xTemp[j])
                        end
                    end
                    break
                end
            end

            stats.rrTime += time() - rrStartTime
                
            if stats.solutionFound                                                                                                                                                                                                
                obj = sum(xTemp[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:nvars)                                                                                                                        
                iterTime = time() - iterStartTime                                                                                                                                                                                 
                if !DEBUG_VERBOSE
                    printRow!(pumpDisplay, stats.pumpIterations, obj, NaN, NaN, 0, 0, iterTime, "randFeasCheck")
                end
                break                                                                                                                                                                                                             
            end
        end

        # if stagnationCount >= DEF_MAX_STAGNATION
        #     stats.restartCount += 1

        #     if stats.restartCount >= DEF_MAX_RESTARTS
        #         stats.exitReason = :restart_limit
        #         break
        #     end

        #     if DEBUG_VERBOSE
        #         println("Cycle detected at iteration $(stats.pumpIterations) (restart #$(stats.restartCount))")
        #         println("Perturbing:")
        #     end

        #     Random.seed!(DEF_RANDOM_SEED + stats.restartCount)  # change seed each restart for reproducibility
        #     perturbSolution!(x, xRound, binIdx, gIntIdx, lpCols)

        #     stagnationCount = 0
        #     bestIntGap = Inf
        #     restarted = true
        # else
        #     xRound .= x  # initialize xRound from LP solution so continuous variables are not left at zero
        #     roundSolution!(xRound, x, intIdx, heur.config.randRound)

        #     if heur.config.lineSearch == :unitary
        #     h = hashSolution(xRound, intIdx)
        #     if h in visitedRounded
        #         perturbSolution!(x, xRound, binIdx, gIntIdx, lpCols)
        #     else
        #         push!(visitedRounded, h)
        #     end
        # end

        xRound .= x  # initialize xRound from LP solution so continuous variables are not left at zero
        roundSolution!(xRound, x, intIdx, heur.config.randRound)
        
        if DEBUG_VERBOSE
            println("Rounding:")
            for i in intIdx
                @printf("  x[%d]: %.3f -> %d\n", i, x[i], Int(xRound[i]))
            end
        end
        
        if heur.config.lineSearch == :unitary
            h = hashSolution(xRound, intIdx)
            if h in visitedRounded
                restarted = true
                stats.restartCount += 1

                if stats.restartCount >= DEF_MAX_RESTARTS
                    stats.exitReason = :restart_limit
                    break
                end

                Random.seed!(DEF_RANDOM_SEED + stats.restartCount)
                perturbSolution!(x, xRound, binIdx, gIntIdx, lpCols)
                h = hashSolution(xRound, intIdx)
            end
            push!(visitedRounded, h)
        end

        # Check if rounded solution is feasible
        if submitSolution(scip, heur_ptr, lpCols, xRound, nvars)
            stats.solutionFound = true
            stats.exitReason = :solutionFound
            result = SCIP.SCIP_FOUNDSOL
            obj = sum(xRound[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:nvars)
            iterTime = time() - iterStartTime
            if !DEBUG_VERBOSE
                printRow!(pumpDisplay, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restart+feasRound" : "feasRound")
            end
            break
        end

        feasible, sol = lpDiving!(scip, lpCols, intIdx, xRound, nvars)                                                                                                                                                              
        if feasible
            if submitSolution(scip, heur_ptr, lpCols, sol, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
                result = SCIP.SCIP_FOUNDSOL
                obj = sum(sol[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:nvars)
                iterTime = time() - iterStartTime
                if !DEBUG_VERBOSE
                    printRow!(pumpDisplay, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restart+divingLP" : "divingLP")
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
            x -> f(x, xRound),
            (storage, x) -> grad!(storage, x, xRound),
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
            x .= xAfter
            xPrev .= xAfter
            continue  # skip rest of checks and go to next FPFW iteration
        else
            xNew = isempty(fwTraj) ? fwResult.x : fwTraj[end][2]
            if heur.config.warmStart && heur.config.fwVariant !== :vanilla
                activeSet = fwResult.active_set
            end
        end

        # Cycle detection: check if we've visited this solution before
        intGap = f(xNew, xRound)

        if isLowerThan(intGap, bestIntGap)
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
                for i in intIdx
                    @printf("  x[%d]: %.3f\n", i, xk[i])
                end
                @printf("Objective = %.3f\n", fk)
            end
        end

        # Step 3: Check feasibility, integrality, distance moved, and objective value
        isIntegral = isSolutionIntegral(xNew, intIdx)
        isFeasible = isSolutionLPFeasible(scip, lpRows, lpCols, xNew, colDict)

        obj = sum(xNew[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:nvars)
        step = f(xNew, xPrev)
        nFrac = count(i -> abs(xNew[i] - round(xNew[i])) > DEF_INT_TOLERANCE, intIdx)
        iterTime = time() - iterStartTime

        if !DEBUG_VERBOSE
            printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, restarted ? "restart" : "")
        end

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isFeasible
            stats.exitReason = :infeasible_fw
            break
        end

        if isIntegral
            # Submit solution to SCIP
            if submitSolution(scip, heur_ptr, lpCols, xNew, nvars)
                stats.solutionFound = true
                stats.exitReason = :solutionFound
                result = SCIP.SCIP_FOUNDSOL

                if DEBUG_VERBOSE
                    println("End solution:")
                    for j in 1:nvars
                        @printf("  x[%d] = %.3f\n", j, xNew[j])
                    end
                end
                break
            else
                # stats.exitReason = :solution_rejected
                # break
                xPrev .= xNew
                x .= xNew
                continue
            end
        else
            # Continue with the projected solution for next FW iteration
            xPrev .= xNew
            x .= xNew
        end
    end

    # Print final summary
    stats.heurTime = time() - startTime
    stats.totalTime = time() - heur.globalStartTime
    stats.primalBound = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    printHeurSummary(stats)

    return (SCIP.SCIP_OKAY, result)
end
