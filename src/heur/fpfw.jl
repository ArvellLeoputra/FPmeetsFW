# Build (or rebuild) the LMO from the current LP and seed its basis
# Used both for the initial build and for the stage-1 -> stage-2 rebuild (freeOld frees the old LPI first)
function setupLMO!(scip, data, config, lpCols, lpRows, colDict, gIntIdx, stage, ncols, nrows; freeOld::Bool=false)
    if config.lmoWarmStart
        if freeOld
            oldLpiRef = Ref(data.lmo.lpi)
            SCIP.@SCIP_CALL SCIP.SCIPlpiFree(oldLpiRef)
        end
        data.lmo = buildLPILMO(scip, lpCols, lpRows, colDict, gIntIdx, config.norm, stage, ncols, nrows, config.verbose)

        newLpiRef = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPgetLPI(scip, newLpiRef)
        if SCIP.SCIPlpiIsOptimal(newLpiRef[]) == SCIP.TRUE
            LPIinitBase(scip, data.lmo, ncols, nrows)
        end
    else
        data.lmo, data.auxConstraintRefs = SCIPbuildLMO(scip, lpCols, lpRows, colDict, gIntIdx, config.norm, ncols, nrows)
    end
end

# Randomized feasibility probe: probabilistically round the current LP point xFrac several times and
# submit each to SCIP. Returns true as soon as one is accepted, leaving the MIP-feasible point in xProbe.
function randFeasCheck!(scip, heur_ptr, lpCols, xFrac, xProbe, intIdx, ncols)
    for _ in 1:DEF_RAND_FEAS_ITER_LIMIT
        xProbe .= xFrac
        for i in intIdx
            # Round up with probability xFrac[i] - floor(xFrac[i]), otherwise round down
            frac = xFrac[i] - floor(xFrac[i])
            xProbe[i] = rand() < frac ? ceil(xFrac[i]) : floor(xFrac[i])
        end

        if submitSolution(scip, heur_ptr, lpCols, xProbe, ncols)
            return true
        end
    end
    return false
end

# Projection using Frank-Wolfe
function fwProject(config, f, grad!, lmo, xFrac, xRound, gIntIdx, intIdx, activeSet, prevGrad, ncols, remainingTime)
    # Rebuilt each call so stateful line searches (e.g. Adaptive) reset per FW solve
    ls = buildLineSearch(config.fwStepSize)

    if config.norm == :manhattan
        # Manhattan start: xFrac plus one aux per general integer, set to |xFrac - xRound| (a feasible start)
        nGInt = length(gIntIdx)
        xStart = zeros(Float64, ncols + nGInt)
        xStart[1:ncols] .= xFrac

        # Set feasible aux values to the distance from xFrac to xRound for each general integer
        for (k, i) in enumerate(gIntIdx)
            xStart[ncols + k] = abs(xFrac[i] - xRound[i])
        end

        gradFn = grad!
    else
        # Smooth norms start: xFrac (a feasible start)
        xStart = xFrac

        # check for grad flips in smooth norms (only when verbose >= 2)
        gradFn = grad!  # default to the original grad! function
        if config.verbose >= 2
            prevGrad .= 0.0
            gradFn = buildGradCheck(grad!, prevGrad, intIdx, xRound)  # wrap grad! to check for flips
        end
    end

    fwResult = runFW(
        config.fwVariant,
        f,
        gradFn,
        lmo;
        x0=xStart,
        activeSet=activeSet,
        warmStart=config.fwWarmStart,
        ls=ls,
        remainingTime=remainingTime,
        fwMaxIterations=config.fwMaxIterations,
        callback=nothing,
        verbose=false
    )

    xProj = config.norm == :manhattan ? fwResult.x[1:ncols] : fwResult.x
    projObj = f(fwResult.x)  # FW objective: distance from xProj to the rounded target
    fwIters = isempty(fwResult.traj_data) ? 0 : fwResult.traj_data[end][1]

    # Update the active set for warm-starting the next FW iteration (if enabled)
    if config.fwWarmStart && config.fwVariant !== :vanilla
        newActiveSet = fwResult.active_set
    else
        newActiveSet = activeSet
    end

    return xProj, projObj, fwIters, newActiveSet
end

# Main FPFW Heuristic Implementation
function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}
    config = heur.config
    data = heur.data
    stats = data.stats

    # Guard the heuristic to only run once per solve
    if data.called > 0
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    # DURINGLPLOOP can fire mid-LP-solve; only proceed once the LP is actually optimal
    if SCIP.SCIPgetLPSolstat(scip) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    # Get LP data
    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    lpCols, lpRows, colDict, binIdx, gIntIdx, initSol = getLPData(scip, ncols, nrows)
    intIdx = [binIdx; gIntIdx]

    # Sanity check: can still legitimately mismatch on an intermediate LP; skip and retry
    dualBound = SCIP.SCIPgetDualbound(scip)
    initObj = origObjective(scip, lpCols, initSol, ncols)
    if SCIP.SCIPisFeasEQ(scip, initObj, dualBound) == SCIP.FALSE
        if config.verbose >= 1
            printstyled("[heuristic skipped]\n", color=:yellow)
            @printf("initObj (%.6f) != dualBound (%.6f) -- LP not final yet, retrying next call\n", initObj, dualBound)
        end
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    data.called += 1

    # Time tracking
    heurStartTime = time()
    scipTime = SCIP.SCIPgetSolvingTime(scip)
    heurTimeLimit = config.timeLimit - scipTime

    # Log initial LP solve info
    if config.verbose >= 1
        lpRootIter = SCIP.SCIPgetNLPIterations(scip)
        nFracVars = SCIP.SCIPgetNLPBranchCands(scip)
        printInitialSolveInfo(initObj, intIdx, scipTime, lpRootIter, nFracVars)
    end

    # Log initial basis info
    if config.verbose >= 2
        cstat, rstat = LPIgetBase(scip, ncols, nrows)
        printInitialBasisInfo(cstat, rstat)
    end

    # Initialize the stage and active integer sets
    # Stage 1: binaries only, Stage 2: all integers
    if countFracVars(scip, binIdx, initSol) > 0
        stage = 1
        activeGIntIdx = Int[]
        activeIntIdx = binIdx
    else
        stage = 2
        activeGIntIdx = gIntIdx
        activeIntIdx = intIdx
    end

    # Per-stage iteration counter (resets at each stage transition)
    stageIter = 0

    if config.verbose >= 2
        printstyled("[debug info]\n", color=:yellow)
    end

    # Build LMO from current LP
    setupLMO!(scip, data, config, lpCols, lpRows, colDict, activeGIntIdx, stage, ncols, nrows)

    # Solution vectors
    xFrac = copy(initSol)              # LP-feasible solution
    prevProj = copy(initSol)           # for distance calculation
    xRound = zeros(Float64, ncols)     # rounded solution (target for FW projection)
    prevRound = zeros(Float64, ncols)  # for cycle detection
    xProbe = zeros(Float64, ncols)     # for randomized feasibility check

    # FW setup
    activeSet = nothing
    prevGrad = zeros(Float64, ncols)  # for detecting gradient flips in smooth norms (only used at verbose >= 2)
    f, grad!, dist = buildFWFunctions(config.norm, binIdx, activeGIntIdx, activeIntIdx, xRound)

    # Perturbation and restart parameters
    avgFlips = max(1, ceil(Int, 0.1 * length(activeIntIdx)))
    visitedRounded = Set{UInt}()  # for cycle detection (only used with unitary step size)
    prevHash = UInt(0)            # for cycle of length 1 detection (only used with unitary step size)

    # Stagnation detection
    bestProjObj = Inf
    stagnationCount = 0
    currentPerturbCount = 0

    # TODO: store the best solution found across iterations, not just the first one
    # foundSolution = nothing

    pumpDisplay = config.verbose == 1 ? setupPumpDisplay() : nothing

    # Main FPFW loop
    result = SCIP.SCIP_DIDNOTFIND
    while true
        # Check time limit
        if timeElapsed(heurStartTime) > heurTimeLimit
            stats.exitReason = TIME_LIMIT
            break
        end

        # Check stage 2 iteration cap
        if stage == 2 && stageIter > DEF_STAGE2_MAX_ITER
            stats.exitReason = ITER_LIMIT
            break
        end

        # Transition to stage 2 if:
        # 1. The solution is binary feasible w.r.t. the binary variables,
        # 2. Stage 1's iteration cap is reached,
        # 3. Stage 1 has already spent DEF_STAGE1_STALL_LIMIT perturb+restart attempts without escaping
        if stage == 1 && (
            countFracVars(scip, binIdx, xFrac) == 0 ||
            stageIter > DEF_STAGE1_MAX_ITER ||
            stats.perturbCount + stats.restartCount >= DEF_STAGE1_STALL_LIMIT
        )
            # Transition to stage 2
            stage = 2
            stageIter = 0
            activeGIntIdx = gIntIdx
            activeIntIdx = intIdx

            # Rebuild the LMO for stage 2 (includes general integers)
            setupLMO!(scip, data, config, lpCols, lpRows, colDict, activeGIntIdx, stage, ncols, nrows; freeOld=true)

            # Reset active set and rebuild the FW functions
            activeSet = nothing
            f, grad!, dist = buildFWFunctions(config.norm, binIdx, activeGIntIdx, activeIntIdx, xRound)

            # Reset perturbation and restart parameters
            avgFlips = max(1, ceil(Int, 0.1 * length(activeIntIdx)))
            empty!(visitedRounded)
            prevHash = UInt(0)

            # Reset stagnation detection
            bestProjObj = Inf
            stagnationCount = 0
            currentPerturbCount = 0
        end

        # Check global iteration limit (should rarely be reached since stage 1 and stage 2 have their own caps)
        if stats.pumpIterations > DEF_MAX_PUMP_ITER
            stats.exitReason = ITER_LIMIT
            break
        end

        stats.pumpIterations += 1
        stageIter += 1
        restarted = false
        perturbed = false
        flips = 0

        iterStartTime = time()  # FP iter start time

        if config.verbose >= 2
            printstyled("[FPFW Iteration $(stats.pumpIterations)]\n"; color=:blue)
        end

        # Random feasibility check (skip the first iteration to save time)
        if config.randFeasCheck && stats.pumpIterations > 1
            rrStartTime = time()
            found = randFeasCheck!(scip, heur_ptr, lpCols, xFrac, xProbe, intIdx, ncols)
            stats.rrTime += timeElapsed(rrStartTime)

            if found
                result = recordSolutionFound!(stats, SOLUTION_RR, scip, heurStartTime)
                origObj = origObjective(scip, lpCols, xProbe, ncols)
                step = dist(xProbe, prevProj)

                logDirectAccept(config, pumpDisplay, stats, stage, origObj, step, heurStartTime, flips, perturbed, restarted, "randFeasCheck", "RandFeasCheck")

                break
            end
        end

        # Step 1: Round LP-feasible solution w.r.t. the current stage's active integer variables
        roundSolution!(xRound, xFrac, activeIntIdx, config.randRound)

        # Rounding debug info
        if config.verbose >= 2
            fracIdx = [i for i in activeIntIdx if !isVarInteger(scip, xFrac[i])]
            nUp = count(i -> xRound[i] > xFrac[i], fracIdx)
            nDown = length(fracIdx) - nUp
            nChanged = count(i -> SCIP.SCIPisEQ(scip, xRound[i], prevRound[i]) == SCIP.FALSE, fracIdx)
            println("  xRound: $nUp up, $nDown down, $nChanged changed / $(length(fracIdx)) fractional")
        end

        # Cycle detection
        if config.fwStepSize == :unitary
            h = hashRounded(xRound, activeIntIdx)
            if h == prevHash
                flips = perturb(scip, xRound, xFrac, binIdx, activeIntIdx, avgFlips, config.verbose >= 2)
                perturbed = flips > 0
                if perturbed
                    stats.perturbCount += 1
                    h = hashRounded(xRound, activeIntIdx)  # rehash after perturbation
                end

            elseif h in visitedRounded
                flips = restart(scip, xRound, xFrac, prevRound, binIdx, activeGIntIdx, lpCols, avgFlips, config.verbose >= 2)
                restarted = flips > 0
                if restarted
                    stats.restartCount += 1
                    h = hashRounded(xRound, activeIntIdx)  # rehash after restart
                    empty!(visitedRounded)  # clear the visited set after a restart
                end
            end

            prevHash = h
            prevRound .= xRound
            push!(visitedRounded, h)
        else
            # Non-unitary FW converges gradually, so exact-hash cycle detection would over-trigger
            # Use objective stagnation instead (no improvement for DEF_MAX_STAGNATION iterations)
            if stagnationCount >= DEF_MAX_STAGNATION
                # currentPerturbCount: total times perturb is called this stage
                if currentPerturbCount < DEF_MAX_PERTURBS
                    currentPerturbCount += 1
                    flips = perturb(scip, xRound, xFrac, binIdx, activeIntIdx, avgFlips, config.verbose >= 2)
                    perturbed = flips > 0
                    if perturbed
                        stats.perturbCount += 1
                        # Reset counters
                        stagnationCount = 0
                        bestProjObj = Inf
                    end

                else  # escalate to a restart once currentPerturbCount reaches DEF_MAX_PERTURBS
                    flips = restart(scip, xRound, xFrac, prevRound, binIdx, activeGIntIdx, lpCols, avgFlips, config.verbose >= 2)
                    restarted = flips > 0
                    if restarted
                        stats.restartCount += 1
                        # Reset counters
                        currentPerturbCount = 0
                        stagnationCount = 0
                        bestProjObj = Inf
                    end
                end
            end
            prevRound .= xRound
        end

        # LMO's rounding-target constraints must reflect the final xRound (post perturb/restart)
        if config.norm == :manhattan
            if config.lmoWarmStart
                LPIupdateRounding!(data.lmo, activeGIntIdx, xRound)
            else
                MOIupdateRounding!(data.lmo, data.auxConstraintRefs, activeGIntIdx, xRound)
            end
        end

        if perturbed && config.verbose >= 3
            println("  xRound = $(xRound[activeIntIdx])  (after perturb)")
        end

        if restarted && config.verbose >= 3
            println("  xRound = $(xRound[activeIntIdx])  (after restart)")
        end

        # Try to submit the rounded solution to SCIP
        if submitSolution(scip, heur_ptr, lpCols, xRound, ncols)
            result = recordSolutionFound!(stats, SOLUTION_ROUND, scip, heurStartTime)
            origObj = origObjective(scip, lpCols, xRound, ncols)
            step = dist(xRound, prevProj)

            logDirectAccept(config, pumpDisplay, stats, stage, origObj, step, heurStartTime, flips, perturbed, restarted, "feasRound", "FeasRound")

            break
        end

        # Skip diveSolve in stage 1 (only binaries are active, so diveSolve can't produce a complete MIP solution)
        if config.useDive && stage == 2
            printstyled("  [diveSolve] fixing integers to the current rounding and diving to solve for continuous values...\n", color=:yellow)
            feasible, sol = diveSolve(scip, lpCols, intIdx, xRound, ncols)
            if feasible
                if submitSolution(scip, heur_ptr, lpCols, sol, ncols)
                    result = recordSolutionFound!(stats, SOLUTION_DIVE, scip, heurStartTime)
                    origObj = origObjective(scip, lpCols, sol, ncols)
                    step = dist(sol, prevProj)

                    logDirectAccept(config, pumpDisplay, stats, stage, origObj, step, heurStartTime, flips, perturbed, restarted, "diveSolve", "DiveSolve")
                    break
                end
            end
        end

        # Step 2: "Projection" using Frank-Wolfe
        remainingTime = heurTimeLimit - timeElapsed(heurStartTime)
        fwStartTime = time()

        xProj, projObj, fwIters, activeSet = fwProject(
            config, f, grad!, data.lmo, xFrac, xRound, activeGIntIdx, activeIntIdx,
            activeSet, prevGrad, ncols, remainingTime
        )

        # Update FW stats
        stats.fwTime += timeElapsed(fwStartTime)
        stats.fwIterations += fwIters

        if config.verbose >= 3
            println("   xProj = $(xProj[activeIntIdx])")
        end

        # Compute metrics for logging
        origObj = origObjective(scip, lpCols, xProj, ncols)
        nFrac = countFracVars(scip, activeIntIdx, xProj)
        step = dist(xProj, prevProj)
        iterTime = timeElapsed(iterStartTime)

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isSolutionLPFeasible(scip, lpRows, lpCols, xProj, colDict)
            stats.exitReason = INFEASIBLE_FW
            logIteration(config, pumpDisplay, stats, stage, origObj, projObj, step, nFrac, fwIters,
                iterTime, heurStartTime, flips, perturbed, restarted, "infeasibleFW", "infeasibleFW")
            break
        end

        # Stagnation tracking (only for non-unitary projections)
        if config.fwStepSize != :unitary
            if SCIP.SCIPisLT(scip, projObj, bestProjObj) == SCIP.TRUE
                if projObj / bestProjObj < 1 - DEF_MIN_IMPROVEMENT
                    stagnationCount = 0
                end
                bestProjObj = projObj
            else
                stagnationCount += 1
            end
        end

        # Step 3: Check feasibility and integrality
        if isSolutionIntegral(scip, xProj, intIdx)
            if submitSolution(scip, heur_ptr, lpCols, xProj, ncols)
                result = recordSolutionFound!(stats, SOLUTION_FWPROJ, scip, heurStartTime)
                logIteration(config, pumpDisplay, stats, stage, origObj, projObj, step, nFrac, fwIters,
                    iterTime, heurStartTime, flips, perturbed, restarted, "feasFWProj", "accepted")

                break
            else  # if SCIP rejects the solution, continue
                logIteration(config, pumpDisplay, stats, stage, origObj, projObj, step, nFrac, fwIters,
                    iterTime, heurStartTime, flips, perturbed, restarted, "rejected", "rejected")
            end
        else  # if xProj is not integral, continue to the next iteration
            logIteration(config, pumpDisplay, stats, stage, origObj, projObj, step, nFrac, fwIters,
                iterTime, heurStartTime, flips, perturbed, restarted, "", "continuing")
        end

        # Continue with the projected solution for next FW iteration
        prevProj .= xProj
        xFrac .= xProj
    end

    # Finalize stats
    stats.heurTime = timeElapsed(heurStartTime)
    stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))
    stats.dualBound = Float64(SCIP.SCIPgetDualbound(scip))
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    # Free the LMO if it was warm-started
    if config.lmoWarmStart && data.lmo !== nothing
        lpiRef = Ref(data.lmo.lpi)
        SCIP.@SCIP_CALL SCIP.SCIPlpiFree(lpiRef)
        data.lmo = nothing
    end

    return (SCIP.SCIP_OKAY, result)
end