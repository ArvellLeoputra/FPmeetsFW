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
    if config.verbose
        cstat, rstat = LPIgetBase(scip, ncols, nrows)
        @printf("Initial basis: cstat L=%d B=%d U=%d | rstat L=%d B=%d U=%d\n",
            count(==(0), cstat), count(==(1), cstat), count(==(2), cstat),
            count(==(0), rstat), count(==(1), rstat), count(==(2), rstat))
    end

    # Build LMO from current LP
    if data.lmo === nothing
        if config.lmoWarmStart
            data.lmo = buildLPILMO(scip, lpCols, lpRows, colDict, intIdx, gIntIdx, config.norm, ncols, nrows, config.verbose)
            LPIinitBase(scip, data.lmo, ncols, nrows)
        else
            data.lmo, data.dConstraintRefs = SCIPbuildLMO(scip, lpCols, lpRows, colDict, gIntIdx, config.norm, ncols, nrows)
        end
    end

    if config.verbose
        printstyled("[debug info]\n", color=:yellow)
    end

    # Solution vectors
    x = copy(initSol)
    xPrev = copy(initSol)  # for distance calculation
    xRound = zeros(Float64, ncols)
    prevRound = zeros(Float64, ncols)
    xTemp = zeros(Float64, ncols)  # for randomized rounding
    
    # FW setup
    activeSet = nothing
    prevGrad = zeros(Float64, ncols)
    f, grad!, dist = buildFWFunctions(config.norm, binIdx, gIntIdx, xRound)
    
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
    
    # Randomized feasibility check parameters
    attempts = min(DEF_RAND_FEAS_ITER_LIMIT, length(intIdx))

    # TODO: store the best solution found across iterations, not just the first one
    # foundSolution = nothing
    
    pumpDisplay = PumpDisplay(PumpDisplayColumn[])
    addColumn!(pumpDisplay, "iter", 6)
    addColumn!(pumpDisplay, "origObj", 15, 2)
    addColumn!(pumpDisplay, "projObj", 15, 4)
    addColumn!(pumpDisplay, "step", 15, 4)
    addColumn!(pumpDisplay, "nFrac", 8)
    addColumn!(pumpDisplay, "fwIters", 10)
    addColumn!(pumpDisplay, "time", 10, 2)
    addColumn!(pumpDisplay, "P", 3)
    addColumn!(pumpDisplay, "R", 3)
    addColumn!(pumpDisplay, "status", 14)

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
            printstyled("[FPFW Iteration $(stats.pumpIterations)]\n"; color=:blue)
        end

        iterStartTime = time()  # FP iter start time

        # Check time limit
        if timeElapsed(heurStartTime) > heurTimeLimit
            stats.exitReason = :time_limit
            break
        end

        # Random feasibility check
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
                origObj = origObjective(scip, lpCols, xTemp, ncols)
                step = dist(xTemp, xPrev)
                elapsed = timeElapsed(heurStartTime)

                if config.verbose
                    @printf("RandFeasCheck: origObj=%.4f step=%.4f elapsed=%.4f\n", origObj, step, elapsed)
                else
                    printRow!(pumpDisplay, stats.pumpIterations, origObj, 0.0, step, 0, 0, elapsed, "", "", "randFeasCheck")
                end

                break
            end
        end

        # Step 1: Round LP feasible solution
        xRound .= x  # initialize xRound from LP solution so continuous variables are not left at zero
        roundSolution!(xRound, x, intIdx, config.randRound)

        if config.norm == :manhattan
            if config.lmoWarmStart
                LPIupdateRounding!(data.lmo, gIntIdx, xRound)
            else
                MOIupdateRounding!(data.lmo, data.dConstraintRefs, gIntIdx, xRound)
            end
        end

        if config.enablePlot
            push!(xIterates, copy(x))
            push!(xRoundIterates, copy(xRound))
        end

        # Cycle detection
        if config.lineSearch == :unitary
            h = hashRounded(xRound, intIdx)
            if h == prevHash
                perturbed = perturb(scip, xRound, x, binIdx, intIdx, avgFlips, config.verbose)
                if perturbed
                    stats.perturbCount += 1
                end

            elseif h in visitedRounded
                restarted = true
                stats.restartCount += 1
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
                    stagnationPerturbCount += 1
                    perturbed = perturb(scip, xRound, x, binIdx, intIdx, avgFlips, config.verbose)
                    if perturbed
                        stats.perturbCount += 1
                    end

                else
                    restarted = true
                    stats.restartCount += 1
                    stagnationPerturbCount = 0
                    restart(scip, xRound, x, prevRound, binIdx, gIntIdx, lpCols, avgFlips)
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

            origObj = origObjective(scip, lpCols, xRound, ncols)
            step = dist(xRound, xPrev)
            elapsed = timeElapsed(heurStartTime)

            if config.verbose
                @printf("FeasRound: origObj=%.4f step=%.4f elapsed=%.4fs P=%s R=%s\n",
                    origObj, step, elapsed, perturbed ? "*" : " ", restarted ? "*" : " ")
            else
                printRow!(pumpDisplay, stats.pumpIterations, origObj, 0.0, step, 0, 0, elapsed, perturbed ? "*" : "", restarted ? "*" : "", "feasRound")
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

                    origObj = origObjective(scip, lpCols, sol, ncols)
                    step = dist(sol, xPrev)
                    elapsed = timeElapsed(heurStartTime)

                    if config.verbose
                        @printf("SubMIPsolve: origObj=%.4f step=%.4f elapsed=%.4fs P=%s R=%s\n",
                            origObj, step, elapsed, perturbed ? "*" : " ", restarted ? "*" : " ")
                    else
                        printRow!(pumpDisplay, stats.pumpIterations, origObj, 0.0, step, 0, 0, elapsed, perturbed ? "*" : "", restarted ? "*" : "", "subMIPsolve")
                    end
                    break
                end
            end
        end

        # Step 2: "Projection" using Frank-Wolfe
        remainingTime = heurTimeLimit - (timeElapsed(heurStartTime))
        fwStartTime = time()
        ls = buildLineSearch(config.lineSearch)

        if config.norm == :manhattan
            # Pad x with d, seeded at the true distance, for a feasible start
            nInt = length(gIntIdx)
            xExt0 = zeros(Float64, ncols + nInt)
            xExt0[1:ncols] .= x

            for (k, i) in enumerate(gIntIdx)
                xExt0[ncols + k] = abs(x[i] - xRound[i])
            end

            fwResult = runFW(
                config.fwVariant,
                f,
                grad!,
                data.lmo;
                x0=xExt0,
                activeSet=nothing,
                warmStart=false,
                ls=ls,
                remainingTime=remainingTime,
                fwMaxIterations=config.fwMaxIterations,
                callback=nothing,
                verbose=false
            )

            intGap = f(fwResult.x)  # sum(d) on the full [x; d] result
            xProj = fwResult.x[1:ncols]
        else
            prevGrad .= 0.0
            # Gradient check function to detect flips in the gradient
            gradCheck! = (storage, x) -> begin
                grad!(storage, x)
                if config.verbose
                    for i in intIdx
                        if prevGrad[i] != 0.0 && storage[i] != prevGrad[i]
                            @printf("  [grad flip] var=%d x=%.6f xRound=%.6f old=%.0f new=%.0f\n",
                                i, x[i], xRound[i], prevGrad[i], storage[i])
                        end
                    end
                end
                copyto!(prevGrad, storage)
                return storage
            end

            fwResult = runFW(
                config.fwVariant,
                f,
                gradCheck!,
                data.lmo;
                x0=x,
                activeSet=activeSet,
                warmStart=config.fwWarmStart,
                ls=ls,
                remainingTime=remainingTime,
                fwMaxIterations=config.fwMaxIterations,
                callback=nothing,
                # Always false: FrankWolfe.jl's own per-iteration verbose output would be very
                # noisy nested inside the outer pump loop. gradCheck! above already provides
                # our own targeted diagnostic (gradient flips) when config.verbose is set.
                verbose=false
            )

            xProj = fwResult.x
            intGap = f(xProj)  # distance to target rounded point (= projObj)
            if config.fwWarmStart && config.fwVariant !== :vanilla
                activeSet = fwResult.active_set
            end
        end

        # Update FW stats
        stats.fwTime += timeElapsed(fwStartTime)
        fwIters = isempty(fwResult.traj_data) ? 0 : fwResult.traj_data[end][1]
        stats.fwIterations += fwIters

        # Compute metrics for logging
        origObj = origObjective(scip, lpCols, xProj, ncols)
        nFrac = countFracVars(scip, intIdx, xProj)
        step = dist(xProj, xPrev)
        iterTime = timeElapsed(iterStartTime)

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isSolutionLPFeasible(scip, lpRows, lpCols, xProj, colDict)
            stats.exitReason = :infeasible_fw
            if config.verbose
                printVerboseIteration(origObj, intGap, step, nFrac, fwIters, iterTime, perturbed, restarted, "infeasibleFW")
            else
                printRow!(pumpDisplay, stats.pumpIterations, origObj, intGap, step, nFrac, fwIters, timeElapsed(heurStartTime), perturbed ? "*" : "", restarted ? "*" : "", "infeasibleFW")
            end
            break
        end

        # Stagnation tracking (only meaningful when the projection is not unitary)
        if config.lineSearch != :unitary
            if SCIP.SCIPisLT(scip, intGap, bestIntGap) == SCIP.TRUE
                bestIntGap = intGap
                stagnationCount = 0
            else
                stagnationCount += 1
            end
        end

        # Step 3: Check feasibility and integrality
        if isSolutionIntegral(scip, xProj, intIdx)
            if submitSolution(scip, heur_ptr, lpCols, xProj, ncols)
                stats.solutionFound = true
                stats.exitReason = :solution_found
                result = SCIP.SCIP_FOUNDSOL
                push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))

                if config.verbose
                    printVerboseIteration(origObj, intGap, step, nFrac, fwIters, iterTime, perturbed, restarted, "accepted")
                else
                    printRow!(pumpDisplay, stats.pumpIterations, origObj, intGap, step, nFrac, fwIters, timeElapsed(heurStartTime), perturbed ? "*" : "", restarted ? "*" : "", "feasFWProj")
                end

                break
            else
                if config.verbose
                    printVerboseIteration(origObj, intGap, step, nFrac, fwIters, iterTime, perturbed, restarted, "rejected")
                else
                    printRow!(pumpDisplay, stats.pumpIterations, origObj, intGap, step, nFrac, fwIters, timeElapsed(heurStartTime), perturbed ? "*" : "", restarted ? "*" : "", "rejected")
                end
            end
        else
            if config.verbose
                printVerboseIteration(origObj, intGap, step, nFrac, fwIters, iterTime, perturbed, restarted, "continuing")
            else
                printRow!(pumpDisplay, stats.pumpIterations, origObj, intGap, step, nFrac, fwIters, timeElapsed(heurStartTime), perturbed ? "*" : "", restarted ? "*" : "", "")
            end
        end
        
        # Continue with the projected solution for next FW iteration
        xPrev .= xProj
        x .= xProj
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