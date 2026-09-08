# Build (or rebuild) the LMO from the current LP and seed its basis
# Used both for the initial build and for the stage-1 -> stage-2 rebuild (freeOld frees the old LPI first)
function setupLMO!(
    scip::Ptr{SCIP.SCIP_},
    data::FPFWRunData,
    config::FPFWConfig,
    lp::LPInfo,
    activeGIntIdx::Vector{Int},
    stage::Int;
    freeOld::Bool = false
)
    if config.lmoWarmStart
        if freeOld
            oldLpiRef = Ref(data.lmo.lpi)
            SCIP.@SCIP_CALL SCIP.SCIPlpiFree(oldLpiRef)
        end

        data.lmo = buildLPILMO(scip, lp.lpCols, lp.lpRows, lp.colDict, activeGIntIdx, config.norm, stage, lp.ncols, lp.nrows, config.verbose)
        newLpiRef = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
        SCIP.@SCIP_CALL SCIP.SCIPgetLPI(scip, newLpiRef)

        if SCIP.SCIPlpiIsOptimal(newLpiRef[]) == SCIP.TRUE
            LPIinitBase(scip, data.lmo, lp.ncols, lp.nrows)
        end
    else
        data.lmo, data.auxConstraintRefs = SCIPbuildLMO(scip, lp.lpCols, lp.lpRows, lp.colDict, activeGIntIdx, config.norm, lp.ncols, lp.nrows)
    end
end

# Check if stage 1 is complete
#  1. All binaries are integral in xFrac, or
#  2. Stage 1 hit its iteration cap, or
#  3. Stage 1 stalled (DEF_STAGE1_NOIMPR_LIMIT iterations with no real bestProjObj gain)
function stage1Complete(
    st::StageState,
    scip::Ptr{SCIP.SCIP_},
    lp::LPInfo,
    xFrac::Vector{Float64}
)::Bool
    return countFracVars(scip, lp.binIdx, xFrac) == 0 ||
           st.stageIter > DEF_STAGE1_MAX_ITER ||
           st.stage1NoImpr > DEF_STAGE1_NOIMPR_LIMIT
end

# Advance `st` from stage 1 (binaries only) to stage 2 (all integers): widen the active sets,
# rebuild the LMO and FW functions, reset the per-stage counters and cycle cache, and resume
# xFrac/prevProj from the closest point stage 1 reached.
function transitionToStage2!(
    st::StageState,
    scip::Ptr{SCIP.SCIP_},
    config::FPFWConfig,
    data::FPFWRunData,
    lp::LPInfo,
    xRound::Vector{Float64},
    xFrac::Vector{Float64},
    prevProj::Vector{Float64},
    visitedRounded::Set{UInt}
)
    st.stage = 2
    st.stageIter = 0
    st.activeGIntIdx = lp.gIntIdx
    st.activeIntIdx = lp.intIdx
    st.avgFlips = max(1, ceil(Int, 0.1 * length(lp.intIdx)))

    setupLMO!(scip, data, config, lp, st.activeGIntIdx, st.stage; freeOld=true)
    f, grad!, dist = buildFWFunctions(config.norm, lp.binIdx, st.activeGIntIdx, st.activeIntIdx, xRound)

    empty!(visitedRounded)
    st.prevHash = UInt(0)
    st.bestProjObj = Inf
    st.stagnationCount = 0
    st.consecutivePerturbs = 0

    # Resume stage 2 from the closest point stage 1 found
    xFrac .= st.closestFrac
    prevProj .= st.closestFrac
    st.closestDist = Inf

    return nothing, f, grad!, dist
end

# Randomized feasibility probe: probabilistically round the current LP point xFrac several times and
# submit each to SCIP. Returns true as soon as one is accepted, leaving the MIP-feasible point in xProbe.
# Called at both stages with all integer constrained variables considered.
function randFeasCheck!(
    scip::Ptr{SCIP.SCIP_},
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
    lp::LPInfo,
    xFrac::Vector{Float64},
    xProbe::Vector{Float64}
)::Bool
    for _ in 1:DEF_RAND_FEAS_ITER_LIMIT
        xProbe .= xFrac
        for i in lp.intIdx
            # Round up with probability xFrac[i] - floor(xFrac[i]), otherwise round down
            frac = xFrac[i] - floor(xFrac[i])
            xProbe[i] = rand() < frac ? ceil(xFrac[i]) : floor(xFrac[i])
        end

        if submitSolution(scip, heur_ptr, lp.lpCols, xProbe, lp.ncols)
            return true
        end
    end
    return false
end

# Projection using Frank-Wolfe
function fwProject(
    config::FPFWConfig,
    f::Function,
    grad!::Function,
    dist::Function,
    lmo::FrankWolfe.LinearMinimizationOracle,
    xFrac::Vector{Float64},
    xRound::Vector{Float64},
    activeGIntIdx::Vector{Int},
    activeIntIdx::Vector{Int},
    activeSet::Union{Nothing, FrankWolfe.ActiveSet},
    prevGrad::Vector{Float64},
    ncols::Int32,
    remainingTime::Float64
)
    # Rebuilt each call so stateful line searches (e.g. Adaptive) reset per FW solve
    ls = buildLineSearch(config.fwStepSize)

    if config.norm == :manhattan
        # Manhattan start: xFrac plus one aux per general integer, set to |xFrac - xRound| (a feasible start)
        nGInt = length(activeGIntIdx)
        xStart = zeros(Float64, ncols + nGInt)
        xStart[1:ncols] .= xFrac

        # Set feasible aux values to the distance from xFrac to xRound for each general integer
        for (k, i) in enumerate(activeGIntIdx)
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
            gradFn = buildGradCheck(grad!, prevGrad, activeIntIdx, xRound)  # wrap grad! to check for flips
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
    projObj = dist(xProj, xRound)  # pure distance to the rounding target

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
    lp = getLPInfo(scip)
    (; lpCols, lpRows, colDict, binIdx, gIntIdx, intIdx, ncols, nrows, initSol) = lp

    data.called += 1

    # Time tracking
    heurStartTime = time()
    rootTime = SCIP.SCIPgetSolvingTime(scip)
    heurTimeLimit = config.timeLimit
    stats.rootTime = rootTime

    # Log initial LP solve info
    if config.verbose >= 1
        initObj = origObjective(scip, lpCols, initSol, ncols)
        printInitialSolveInfo(scip, initObj, intIdx)
    end

    # Log initial basis info
    if config.verbose >= 2
        cstat, rstat = LPIgetBase(scip, ncols, nrows)
        printInitialBasisInfo(cstat, rstat)
    end

    # Per-stage state
    if countFracVars(scip, binIdx, initSol) > 0
        st = StageState(; stage=1, activeIntIdx=binIdx, activeGIntIdx=Int[], closestFrac=copy(initSol))
    else
        st = StageState(; stage=2, activeIntIdx=intIdx, activeGIntIdx=gIntIdx, closestFrac=copy(initSol))
    end

    if config.verbose >= 2
        printstyled("[debug info]\n", color=:yellow)
    end

    # Build LMO from current LP
    setupLMO!(scip, data, config, lp, st.activeGIntIdx, st.stage)

    # Solution vectors
    xFrac = copy(initSol)              # LP-feasible solution
    prevProj = copy(initSol)           # for distance calculation
    xRound = zeros(Float64, ncols)     # rounded solution (target for FW projection)
    prevRound = zeros(Float64, ncols)  # for cycle detection
    xProbe = zeros(Float64, ncols)     # for randomized feasibility check

    # FW setup
    activeSet = nothing
    prevGrad = zeros(Float64, ncols)  # for detecting gradient flips in smooth norms (only used at verbose >= 2)
    f, grad!, dist = buildFWFunctions(config.norm, binIdx, st.activeGIntIdx, st.activeIntIdx, xRound)

    # Objective feasibility pump: initial weight, disabled when there's no objective
    alpha = lp.objScale > 0.0 ? config.alpha : 0.0

    # Cycle detection cache (only used with unitary step size)
    visitedRounded = Set{UInt}()

    # TODO: store the best solution found across iterations, not just the first one
    # foundSolution = nothing

    pumpDisplay = config.verbose == 1 ? setupPumpDisplay(config) : nothing

    # Main FPFW loop
    result = SCIP.SCIP_DIDNOTFIND
    while true
        # Check time limit
        if timeElapsed(heurStartTime) > heurTimeLimit
            stats.exitReason = TIME_LIMIT
            break
        end

        # Check stage 2 iteration cap
        if st.stage == 2 && st.stageIter > DEF_STAGE2_MAX_ITER
            stats.exitReason = ITER_LIMIT
            break
        end

        # Advance to stage 2 once stage 1 is done (binaries satisfied / iter cap / stalled)
        if st.stage == 1 && stage1Complete(st, scip, lp, xFrac)
            activeSet, f, grad!, dist = transitionToStage2!(st, scip, config, data, lp,
                                                            xRound, xFrac, prevProj, visitedRounded)
        end

        # Check global iteration limit (should rarely be reached since stage 1 and stage 2 have their own caps)
        if stats.pumpIterations > DEF_MAX_PUMP_ITER
            stats.exitReason = ITER_LIMIT
            break
        end

        stats.pumpIterations += 1
        st.stageIter += 1
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
            found = randFeasCheck!(scip, heur_ptr, lp, xFrac, xProbe)
            stats.rrTime += timeElapsed(rrStartTime)

            if found
                result = recordSolutionFound!(stats, SOLUTION_RR, scip, heurStartTime)
                origObj = origObjective(scip, lpCols, xProbe, ncols)
                step = dist(xProbe, prevProj)

                logDirectAccept(config, pumpDisplay, stats, st.stage, alpha, origObj, step, heurStartTime, flips, perturbed, restarted, "randFeasCheck", "RandFeasCheck")

                break
            end
        end

        # Step 1: Round LP-feasible solution w.r.t. the current stage's active integer variables
        roundSolution!(xRound, xFrac, st.activeIntIdx, config.randRound)

        # Rounding debug info
        if config.verbose >= 2
            fracIdx = [i for i in st.activeIntIdx if !isVarInteger(scip, xFrac[i])]
            nUp = count(i -> xRound[i] > xFrac[i], fracIdx)
            nDown = length(fracIdx) - nUp
            nChanged = count(i -> SCIP.SCIPisEQ(scip, xRound[i], prevRound[i]) == SCIP.FALSE, fracIdx)
            println("  xRound: $nUp up, $nDown down, $nChanged changed / $(length(fracIdx)) fractional")
        end

        # Cycle / stagnation detection
        if config.fwStepSize == :unitary
            h = hashRounded(xRound, st.activeIntIdx)

            # stucked: the rounded solution is identical to the previous iteration's
            # cycled: the rounded solution has been seen before (but not in the previous iteration)
            stucked = h == st.prevHash
            cycled = !stucked && h in visitedRounded
            stagnated = st.stagnationCount >= DEF_MAX_STAGNATION

            if stucked || cycled || stagnated
                # Restart (rather than perturb) on a genuine longer cycle or if the perturbation limit has been reached
                doRestart = cycled || st.consecutivePerturbs >= DEF_MAX_PERTURBS

                if !doRestart
                    flips = perturb(scip, xRound, xFrac, binIdx, st.activeIntIdx, st.avgFlips, config.verbose >= 2)
                    perturbed = flips > 0
                    if perturbed
                        stats.perturbCount += 1
                        st.consecutivePerturbs += 1
                        st.stagnationCount = 0
                        st.bestProjObj = Inf
                        h = hashRounded(xRound, st.activeIntIdx)  # rehash after perturbation
                    else
                        # edge case: perturbation failed to flip any variables, so escalate to a restart
                        # happens only when the LP solution is already integral but rejected by SCIP
                        doRestart = true
                    end
                end

                if doRestart
                    flips = restart(scip, xRound, xFrac, prevRound, binIdx, st.activeGIntIdx, lpCols, st.avgFlips, config.verbose >= 2)
                    restarted = flips > 0
                    if restarted
                        stats.restartCount += 1
                        st.consecutivePerturbs = 0
                        st.stagnationCount = 0
                        st.bestProjObj = Inf
                        h = hashRounded(xRound, st.activeIntIdx)  # rehash after restart
                        empty!(visitedRounded)  # clear the visited set after a restart
                    end
                end
            end

            st.prevHash = h
            prevRound .= xRound
            push!(visitedRounded, h)
        else
            # Non-unitary FW converges gradually, so exact-hash cycle detection would over-trigger
            # Use objective stagnation instead (no improvement for DEF_MAX_STAGNATION iterations)
            if st.stagnationCount >= DEF_MAX_STAGNATION
                # consecutivePerturbs: perturbs since the last restart (or stage start)
                if st.consecutivePerturbs < DEF_MAX_PERTURBS
                    st.consecutivePerturbs += 1
                    flips = perturb(scip, xRound, xFrac, binIdx, st.activeIntIdx, st.avgFlips, config.verbose >= 2)
                    perturbed = flips > 0
                    if perturbed
                        stats.perturbCount += 1
                        # Reset counters
                        st.stagnationCount = 0
                        st.bestProjObj = Inf
                    end

                else  # escalate to a restart once consecutivePerturbs reaches DEF_MAX_PERTURBS
                    flips = restart(scip, xRound, xFrac, prevRound, binIdx, st.activeGIntIdx, lpCols, st.avgFlips, config.verbose >= 2)
                    restarted = flips > 0
                    if restarted
                        stats.restartCount += 1
                        # Reset counters
                        st.consecutivePerturbs = 0
                        st.stagnationCount = 0
                        st.bestProjObj = Inf
                    end
                end
            end
            prevRound .= xRound
        end

        # LMO's rounding-target constraints must reflect the final xRound (post perturb/restart)
        if config.norm == :manhattan
            if config.lmoWarmStart
                LPIupdateRounding!(data.lmo, st.activeGIntIdx, xRound)
            else
                MOIupdateRounding!(data.lmo, data.auxConstraintRefs, st.activeGIntIdx, xRound)
            end
        end

        if perturbed && config.verbose >= 4
            println("  xRound = $(xRound[st.activeIntIdx])  (after perturb)")
        end

        if restarted && config.verbose >= 4
            println("  xRound = $(xRound[st.activeIntIdx])  (after restart)")
        end

        # Try to submit the rounded solution to SCIP
        if submitSolution(scip, heur_ptr, lpCols, xRound, ncols)
            result = recordSolutionFound!(stats, SOLUTION_ROUND, scip, heurStartTime)
            origObj = origObjective(scip, lpCols, xRound, ncols)
            step = dist(xRound, prevProj)

            logDirectAccept(config, pumpDisplay, stats, st.stage, alpha, origObj, step, heurStartTime, flips, perturbed, restarted, "feasRound", "FeasRound")

            break
        end

        # Skip diveSolve in stage 1 (only binaries are active, so diveSolve can't produce a complete MIP solution)
        if config.useDive && st.stage == 2
            printstyled("  [diveSolve] fixing integers to the current rounding and diving to solve for continuous values...\n", color=:yellow)
            feasible, sol = diveSolve(scip, lpCols, intIdx, xRound, ncols)
            if feasible
                if submitSolution(scip, heur_ptr, lpCols, sol, ncols)
                    result = recordSolutionFound!(stats, SOLUTION_DIVE, scip, heurStartTime)
                    origObj = origObjective(scip, lpCols, sol, ncols)
                    step = dist(sol, prevProj)

                    logDirectAccept(config, pumpDisplay, stats, st.stage, alpha, origObj, step, heurStartTime, flips, perturbed, restarted, "diveSolve", "DiveSolve")
                    break
                end
            end
        end

        # Step 2: "Projection" using Frank-Wolfe
        distScale = sqrt(length(st.activeIntIdx))  # || delta ||
        weight = alpha > 0.0 ? distScale / max(lp.objScale, 1.0) : 0.0
        fOFP, gradOFP! = buildOFPFunctions(f, grad!, alpha, weight, lp.objCoeffs, ncols)

        remainingTime = heurTimeLimit - timeElapsed(heurStartTime)
        fwStartTime = time()

        # Bound the LMO's own LP solves to what's left of the budget, so a single slow LP can't run unbounded
        data.lmo.deadline[] = time() + remainingTime

        # Pre-declared: a `try` block inside a `while` loop is its own scope in Julia, so
        # first-assigning these inside the try below would not be visible after it.
        local xProj, projObj, fwIters
        try
            xProj, projObj, fwIters, activeSet = fwProject(
                config, fOFP, gradOFP!, dist, data.lmo, xFrac, xRound, st.activeGIntIdx, st.activeIntIdx,
                activeSet, prevGrad, ncols, remainingTime
            )
        catch e
            e isa LMODeadlineExceeded || rethrow()
            stats.fwTime += timeElapsed(fwStartTime)
            stats.exitReason = TIME_LIMIT
            break
        end

        # Update FW stats
        stats.fwTime += timeElapsed(fwStartTime)
        stats.fwIterations += fwIters

        if config.verbose >= 4
            println("   xProj = $(xProj[st.activeIntIdx])")
        end

        # Compute metrics for logging
        origObj = origObjective(scip, lpCols, xProj, ncols)
        nFrac = countFracVars(scip, st.activeIntIdx, xProj)
        step = dist(xProj, prevProj)
        iterTime = timeElapsed(iterStartTime)

        # Safety check: FW must always return a feasible point (LP polytope is preserved)
        if !isSolutionLPFeasible(scip, lpRows, lpCols, xProj, colDict)
            stats.exitReason = INFEASIBLE_FW
            logIteration(config, pumpDisplay, stats, st.stage, alpha, origObj, projObj, step, nFrac, fwIters,
                iterTime, heurStartTime, flips, perturbed, restarted, "infeasibleFW", "infeasibleFW")
            break
        end

        # Stagnation tracking
        # stagnationCount: perturb/restart
        # stage1NoImpr: stage-1 stall exit
        if SCIP.SCIPisLT(scip, projObj, st.bestProjObj) == SCIP.TRUE
            if projObj / st.bestProjObj < 1 - DEF_MIN_IMPROVEMENT
                st.stagnationCount = 0
                st.stage == 1 && (st.stage1NoImpr = 0)
            end
            st.bestProjObj = projObj
        else
            st.stagnationCount += 1
            st.stage == 1 && (st.stage1NoImpr += 1)
        end

        # Closest point of the stage
        if SCIP.SCIPisLT(scip, projObj, st.closestDist) == SCIP.TRUE
            st.closestDist = projObj
            st.closestFrac .= xProj
        end

        # Step 3: Check feasibility and integrality
        if isSolutionIntegral(scip, xProj, intIdx)
            if submitSolution(scip, heur_ptr, lpCols, xProj, ncols)
                result = recordSolutionFound!(stats, SOLUTION_FWPROJ, scip, heurStartTime)
                logIteration(config, pumpDisplay, stats, st.stage, alpha, origObj, projObj, step, nFrac, fwIters,
                    iterTime, heurStartTime, flips, perturbed, restarted, "feasFWProj", "accepted")

                break
            else  # if SCIP rejects the solution, continue
                logIteration(config, pumpDisplay, stats, st.stage, alpha, origObj, projObj, step, nFrac, fwIters,
                    iterTime, heurStartTime, flips, perturbed, restarted, "rejected", "rejected")
            end
        else  # if xProj is not integral, continue to the next iteration
            logIteration(config, pumpDisplay, stats, st.stage, alpha, origObj, projObj, step, nFrac, fwIters,
                iterTime, heurStartTime, flips, perturbed, restarted, "", "continuing")
        end

        if alpha > 0.0
            alpha *= config.alphaFactor
            if alpha <= DEF_ALPHA_MIN
                alpha = 0.0
            end
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