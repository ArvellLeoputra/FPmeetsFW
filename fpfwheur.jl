# Main FPFW Heuristic Implementation
function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}
    data = heur.data
    stats = data.stats
    config = heur.config

    data.called += 1
    if data.called > 1
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    result = SCIP.SCIP_DIDNOTFIND

    heurStartTime = time()
    scipTime = SCIP.SCIPgetSolvingTime(scip)
    heurTimeLimit = config.timeLimit - scipTime

    # Get LP data
    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    lpCols, lpRows, colDict, binIdx, gIntIdx, initSol = getLPData(scip, ncols, nrows)
    intIdx = [binIdx; gIntIdx]

    # Log initial LP solve info
    stats.dualBound = SCIP.SCIPgetDualbound(scip)
    lpRootIter = SCIP.SCIPgetNLPIterations(scip)
    nFracVars = SCIP.SCIPgetNLPBranchCands(scip)

    printstyled("[initialSolve]\n", color=:cyan)
    @printf("Initial LP: lpiter=%d obj=%.2f frac=%d/%d time=%.2fs\n", 
            lpRootIter, stats.dualBound, nFracVars, length(intIdx), scipTime)

    # Build LMO from current LP
    if data.lmo === nothing
        if config.lmoWarmStart
            data.lmo = buildLPILMO(scip, lpCols, lpRows, colDict, ncols, nrows, config.verbose)
            LPIinitBase(scip, data.lmo)
        else
            data.lmo = SCIPbuildLMO(scip, lpCols, lpRows, colDict, ncols, nrows)
        end
    end

    if config.verbose
        printstyled("[debug info]\n", color=:yellow)
        cstat, rstat = LPIgetBase(scip, ncols, nrows)
        @printf("Initial basis: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d\n",
            count(==(0), cstat), count(==(1), cstat), count(==(2), cstat),
            count(==(0), rstat), count(==(1), rstat), count(==(2), rstat))
    end

    # Solution vectors
    x = copy(initSol)
    xPrev = copy(initSol)  # for distance calculation
    xRound = zeros(Float64, ncols)
    xTemp = zeros(Float64, ncols)  # for randomized rounding
    prevRound = zeros(Float64, ncols)
    # xAfter = zeros(Float64, ncols)  # for post-step x in callback
    # xRoundEscape = zeros(Float64, ncols)

    # FW setup
    activeSet = nothing
    f, grad!, dist = buildFWFunctions(config.norm, intIdx, xRound)
    # fwCallback = (state, args...) -> begin
    #     # Skip FW bookkeeping steps where d and gamma are not meaningful
    #     if state.step_type === FrankWolfe.ST_LAST || state.step_type === FrankWolfe.ST_POSTPROCESS
    #         return true
    #     end
    #
    #     if state.d === nothing || state.gamma === nothing
    #         return true
    #     end
    #
    #     xAfter .= state.x .- state.gamma .* state.d
    #     fwStepCount += 1
    #
    #     # Check if xAfter rounds to a different target than xRound
    #     # If so, stop FW early and use xAfter as the new starting point
    #     if DEF_FW_ESCAPE
    #         roundSolution!(xRoundEscape, xAfter, intIdx, config.randRound)
    #         if !areSolutionsEqual(scip, intIdx, xRoundEscape, xRound)
    #             fwEscaped = true
    #
    #             obj = sum(xAfter[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
    #             step = dist(xAfter, xPrev)
    #             intGap = f(xAfter)
    #             nFrac = countFracVars(scip, intIdx, xAfter)
    #             escFwIters = state.t
    #             iterTime = timeElapsed(iterStartTime)
    #
    #             if config.verbose
    #                 @printf("FW escaped rounding target at step %d: obj=%.4f projObj=%.4f step=%.4f nFrac=%d\n",
    #                         escFwIters, obj, intGap, step, nFrac)
    #             else
    #                 printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, escFwIters, iterTime, "escape")
    #             end
    #
    #             return false  # stop FW iter early
    #         end
    #     end
    #
    #     return true  # continue FW iter
    # end
    
    # Iterate collection for 2D plotting
    xIterates = Vector{Vector{Float64}}()
    xRoundIterates = Vector{Vector{Float64}}()

    # Cycle detection
    restarted = false
    perturbed = false
    avgFlips = max(1, ceil(Int, 0.1 * length(intIdx)))
    visitedRounded = Set{UInt}()
    prevHash = UInt(0)

    # Stagnation detection
    bestIntGap = Inf
    stagnationCount = 0
    stagnationPerturbCount = 0
    
    # FW escape flag
    # fwEscaped = false
    # fwStepCount = 0
    
    # Randomized feasibility check parameters
    attempts = min(DEF_RAND_FEAS_ITER_LIMIT, length(intIdx))

    # TODO: store the best solution found across iterations, not just the first one
    # foundSolution = nothing
    
    pumpDisplay = PumpDisplay(PumpDisplayColumn[])
    addColumn!(pumpDisplay, "pumpIter", 10)
    addColumn!(pumpDisplay, "obj", 15, 4)
    addColumn!(pumpDisplay, "projObj", 15, 4)
    addColumn!(pumpDisplay, "step", 15, 4)
    addColumn!(pumpDisplay, "nFrac", 8)
    addColumn!(pumpDisplay, "fwIters", 10)
    addColumn!(pumpDisplay, "time", 10, 2)
    addColumn!(pumpDisplay, "flag", 18)

    if !config.verbose
        printstyled("[pump]\n", color=:cyan)
        printHeader!(pumpDisplay)
    end

    # Main FPFW loop
    while true
        stats.pumpIterations += 1
        restarted = false
        perturbed = false

        if config.verbose
            printstyled("\nFPFW Iteration $(stats.pumpIterations)\n"; color=:blue)
        end

        iterStartTime = time()  # FP iter start time

        # Check time limit
        if timeElapsed(heurStartTime) > heurTimeLimit
            stats.exitReason = :time_limit
            break
        end

        if config.randFeasCheck && stats.pumpIterations > 1  # skip randomized rounding in the first iteration to save time
            rrStartTime = time()
            for _ in 1:attempts
                xTemp .= x
                for i in intIdx
                    frac = x[i] - floor(x[i])
                    xTemp[i] = rand() < frac ? ceil(x[i]) : floor(x[i])
                end

                if submitSolution(scip, heur_ptr, lpCols, xTemp, ncols)
                    stats.solutionFound = true
                    stats.exitReason = :rr_solution_found
                    result = SCIP.SCIP_FOUNDSOL
                    push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))
                    break
                end
            end

            stats.rrTime += timeElapsed(rrStartTime)

            if stats.solutionFound
                obj = sum(xTemp[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
                iterTime = timeElapsed(iterStartTime)
                if !config.verbose
                    printRow!(pumpDisplay, stats.pumpIterations, obj, NaN, NaN, 0, 0, iterTime, "randFeasCheck")
                end
                break
            end
        end

        xRound .= x  # initialize xRound from LP solution so continuous variables are not left at zero
        roundSolution!(xRound, x, intIdx, config.randRound)
        push!(xIterates, copy(x))
        push!(xRoundIterates, copy(xRound))

        if config.lineSearch == :unitary
            h = hashRounded(xRound, intIdx)
            if h == prevHash
                perturbed = true
                stats.perturbCount += 1
                
                if config.verbose
                    println("Cycle detected at iteration $(stats.pumpIterations) (perturb #$(stats.perturbCount))")
                end
                
                perturb(scip, xRound, x, binIdx, intIdx, avgFlips, config.verbose)

            elseif h in visitedRounded
                restarted = true
                stats.restartCount += 1

                if config.verbose
                    println("Cycle detected at iteration $(stats.pumpIterations) (restart #$(stats.restartCount))")
                end
                
                restart(scip, xRound, x, prevRound, binIdx, gIntIdx, lpCols, avgFlips)
                empty!(visitedRounded)
            end

            h = hashRounded(xRound, intIdx)
            prevHash = h
            prevRound .= xRound
            push!(visitedRounded, h)
        else
            if stagnationCount >= DEF_MAX_STAGNATION
                if stagnationPerturbCount < DEF_STAGNATION_RESTART_THRESHOLD
                    perturbed = true
                    stats.perturbCount += 1
                    stagnationPerturbCount += 1
                    perturb(scip, xRound, x, binIdx, intIdx, avgFlips, config.verbose)

                    if config.verbose
                        println("Stagnation detected at iteration $(stats.pumpIterations) (perturb #$(stats.perturbCount))")
                    end

                else
                    restarted = true
                    stats.restartCount += 1
                    stagnationPerturbCount = 0
                    restart(scip, xRound, x, prevRound, binIdx, gIntIdx, lpCols, avgFlips)

                    if config.verbose
                        println("More stagnations detected at iteration $(stats.pumpIterations) (restart #$(stats.restartCount))")
                    end
                end
                prevRound .= xRound
                stagnationCount = 0
                bestIntGap = Inf
            end
        end

        # Check if rounded solution is feasible
        if submitSolution(scip, heur_ptr, lpCols, xRound, ncols)
            stats.solutionFound = true
            stats.exitReason = :solution_found
            result = SCIP.SCIP_FOUNDSOL
            push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))

            obj = sum(xRound[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
            iterTime = timeElapsed(iterStartTime)
            if !config.verbose
                printRow!(pumpDisplay, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restarted+feasRound" : perturbed ? "perturbed+feasRound" : "feasRound")
            end
            break
        end

        if config.useSubMIP
            feasible, sol = subMIPsolve(scip, lpCols, intIdx, xRound, ncols)
            if feasible
                if submitSolution(scip, heur_ptr, lpCols, sol, ncols)
                    stats.solutionFound = true
                    stats.exitReason = :solution_found
                    result = SCIP.SCIP_FOUNDSOL
                    push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))

                    obj = sum(sol[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
                    iterTime = timeElapsed(iterStartTime)
                    if !config.verbose
                        printRow!(pumpDisplay, stats.pumpIterations, obj, 0.0, NaN, 0, 0, iterTime, restarted ? "restarted+subMIPsolve" : perturbed ? "perturbed+subMIPsolve" : "subMIPsolve")
                    end
                    break
                end
            end
        end

        # fwEscaped = false
        # fwStepCount = 0

        # Step 2: "Projection" using Frank-Wolfe
        remainingTime = heurTimeLimit - (timeElapsed(heurStartTime))
        fwStartTime = time()
        ls = buildLineSearch(config.lineSearch)

        fwResult = run_fw(
            config.fwVariant,
            f,
            grad!,
            data.lmo,
            x,
            activeSet,
            config.fwWarmStart,
            ls,
            nothing,
            remainingTime,
            config.fwMaxIterations,
            false
        )

        stats.fwTime += timeElapsed(fwStartTime)
        fwIters = isempty(fwResult.traj_data) ? 0 : fwResult.traj_data[end][1]
        stats.fwIterations += fwIters

        # If FW gives us intermediate solutions via callback that escape the rounding target, we immediately jump to next iteration with the escaped solution
        # if fwEscaped
        #     x .= xAfter
        #     xPrev .= xAfter
        #     continue  # skip rest of checks and go to next FPFW iteration
        # else
        xProj = fwResult.x
        if config.fwWarmStart && config.fwVariant !== :vanilla
            activeSet = fwResult.active_set
        end
        # end

        # Cycle detection: check if we've visited this solution before
        intGap = f(xProj)

        if config.lineSearch != :unitary
            if SCIP.SCIPisLT(scip, intGap, bestIntGap) == SCIP.TRUE
                bestIntGap = intGap
                stagnationCount = 0
            else
                stagnationCount += 1
            end
        end


        # Step 3: Check feasibility, integrality, distance moved, and objective value
        isIntegral = isSolutionIntegral(scip, xProj, intIdx)
        isFeasible = isSolutionLPFeasible(scip, lpRows, lpCols, xProj, colDict)

        obj = sum(xProj[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lpCols[j])) for j in 1:ncols)
        step = dist(xProj, xPrev)
        nFrac = countFracVars(scip, intIdx, xProj)
        iterTime = timeElapsed(iterStartTime)

        if config.verbose
            @printf("Projection: projObj=%.4f nFrac=%d obj=%.4f\n", intGap, nFrac, obj)
        end

        flag = restarted ? "restarted" : perturbed ? "perturbed" : ""

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isFeasible
            stats.exitReason = :infeasible_fw
            if !config.verbose
                printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, restarted ? "restarted+infeasible_fw" : perturbed ? "perturbed+infeasible_fw" : "infeasible_fw")
            end
            break
        end

        if isIntegral
            # Submit solution to SCIP
            if submitSolution(scip, heur_ptr, lpCols, xProj, ncols)
                stats.solutionFound = true
                stats.exitReason = :solution_found
                result = SCIP.SCIP_FOUNDSOL
                push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))

                if !config.verbose
                    printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, flag)
                end
                break
            else
                flag = restarted ? "restarted+rejected" : perturbed ? "perturbed+rejected" : "rejected"
                if !config.verbose
                    printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, flag)
                end
                xPrev .= xProj
                x .= xProj
                continue
            end
        else
            # Continue with the projected solution for next FW iteration
            if !config.verbose
                printRow!(pumpDisplay, stats.pumpIterations, obj, intGap, step, nFrac, fwIters, iterTime, flag)
            end
            xPrev .= xProj
            x .= xProj
        end
    end

    stats.heurTime = timeElapsed(heurStartTime)
    stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    if config.lmoWarmStart && data.lmo !== nothing
        lpiRef = Ref(data.lmo.lpi)
        SCIP.@SCIP_CALL SCIP.SCIPlpiFree(lpiRef)
    end

    if config.enablePlot && length(intIdx) == 2
        instanceName = unsafe_string(SCIP.SCIPgetProbName(scip))
        plotIterates(scip, lpRows, lpCols, colDict, intIdx, xIterates, xRoundIterates, config, instanceName)
    end

    return (SCIP.SCIP_OKAY, result)
end
