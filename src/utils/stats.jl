function timeElapsed(start::Float64)
    return time() - start
end

function computePrimalIntegral(events::Vector{Tuple{Float64, Float64}}, totalTime::Float64)
    area = 0.0
    prevTime = 0.0
    prevGap  = 1.0

    for (t, gap) in events
        area += prevGap * (t - prevTime)
        prevTime = t
        prevGap = gap
    end

    area += prevGap * (totalTime - prevTime)
    return area
end

function buildStats(scip::SCIP.SCIPData, heur::FPFWHeuristic, totalTime::Float64)
    stats = heur.data.stats
    if heur.data.called > 0
        stats.totalTime = totalTime
        stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)
        return stats
    end

    stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))
    stats.dualBound = Float64(SCIP.SCIPgetDualbound(scip))
    stats.gap = Float64(SCIP.SCIPgetGap(scip))

    stats.totalTime = totalTime
    stats.solutionFound = SCIP.SCIPgetNSols(scip) > 0

    status = SCIP.SCIPgetStatus(scip)
    if status == SCIP.SCIP_STATUS_OPTIMAL
        stats.exitReason = SCIP_OPTIMAL
    elseif status == SCIP.SCIP_STATUS_INFEASIBLE
        stats.exitReason = SCIP_INFEASIBLE
    elseif status == SCIP.SCIP_STATUS_UNBOUNDED
        stats.exitReason = SCIP_UNBOUNDED
    elseif status == SCIP.SCIP_STATUS_TIMELIMIT
        stats.exitReason = SCIP_TIME_LIMIT
    elseif status == SCIP.SCIP_STATUS_NODELIMIT
        stats.exitReason = SCIP_NODE_LIMIT
    else
        stats.exitReason = SCIP_UNKNOWN
    end

    stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)

    return stats
end

function recordSolutionFound!(stats, exitReason, scip, heurStartTime)
    stats.solutionFound = true
    stats.exitReason = exitReason
    push!(stats.primalEvents, (timeElapsed(heurStartTime), min(1.0, Float64(SCIP.SCIPgetGap(scip)))))
    return SCIP.SCIP_FOUNDSOL
end

function printRunInfo(scip::SCIP.SCIPData, config::FPFWConfig)
    name = unsafe_string(SCIP.SCIPgetProbName(scip))
    scipVersion = "$(SCIP.SCIPmajorVersion()).$(SCIP.SCIPminorVersion()).$(SCIP.SCIPtechVersion())"
    varCount = SCIP.SCIPgetNOrigVars(scip)
    binCount = SCIP.SCIPgetNBinVars(scip)
    intCount = SCIP.SCIPgetNIntVars(scip)
    contCount = SCIP.SCIPgetNContVars(scip)

    printstyled("[run info]\n", color=:cyan)
    println("runName = $(config.runName)")
    println("scipVersion = $scipVersion")
    println("instance = $name")
    println("totalVars = $varCount")
    println("binaryVars = $binCount")
    println("integerVars = $intCount")
    println("continuousVars = $contCount")
end

function printConfigs(config::FPFWConfig)
    printstyled("[FPFW configs]\n", color=:cyan)
    println("norm = $(config.norm)")
    println("fwVariant = $(config.fwVariant)")
    println("fwMaxIterations = $(config.fwMaxIterations)")
    println("fwStepSize = $(config.fwStepSize)")
    println("timeLimit = $(config.timeLimit)")
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
    println("fwWarmStart = $(config.fwWarmStart ? "enabled" : "disabled")")
    println("lmoWarmStart = $(config.lmoWarmStart ? "enabled" : "disabled")")
    println("useDive = $(config.useDive ? "enabled" : "disabled")")
    println("presolve = $(config.presolve ? "enabled" : "disabled")")
    println("seed = $(config.seed)")
    println("verbose = $(config.verbose)")
end

function printResults(stats::FPFWStats)
    exitMsg = if stats.exitReason == TIME_LIMIT
        "time limit reached"
    elseif stats.exitReason == ITER_LIMIT
        "iteration limit reached without a solution"
    elseif stats.exitReason == INFEASIBLE_FW
        "FW returned a point outside the feasible polytope (numerical error)"
    elseif stats.exitReason == SOLUTION_ROUND
        "integer feasible solution found by direct rounding at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == SOLUTION_DIVE
        "integer feasible solution found by dive (integers fixed, LP solved) at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == SOLUTION_FWPROJ
        "integer feasible solution found by FW projection at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == SOLUTION_RR
        "integer feasible solution found by randomized rounding at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == SCIP_OPTIMAL
        "problem solved to optimality by SCIP before heuristic was called"
    elseif stats.exitReason == SCIP_INFEASIBLE
        "LP is infeasible, thus MIP is infeasible"
    elseif stats.exitReason == SCIP_TIME_LIMIT
        "time limit reached before heuristic was called"
    elseif stats.exitReason == SCIP_NODE_LIMIT
        "SCIP node limit reached before heuristic was called"
    elseif stats.exitReason == SCIP_UNBOUNDED
        "LP is unbounded, thus MIP is unbounded"
    else
        "unknown exit"
    end

    printstyled("[result]\n", color=:cyan)
    println("primalBound = $(round(stats.primalBound, digits=4))")
    println("dualBound = $(round(stats.dualBound, digits=4))")
    gapStr = stats.gap >= 1e19 ? @sprintf("%.0e", stats.gap) : @sprintf("%.2f %%", stats.gap * 100)
    println("gap = $gapStr")
    println("primalIntegral = $(round(stats.primalIntegral, digits=4))")
    println("solFound = $(stats.solutionFound)")
    println("totalTime = $(round(stats.totalTime, digits=2))s")

    if !startswith(string(stats.exitReason), "SCIP_")
        println("totalHeurTime = $(round(stats.heurTime, digits=2))s")
        println("fwTime = $(round(stats.fwTime, digits=2))s")
        println("randRoundTime = $(round(stats.rrTime, digits=2))s")
        println("pumpIterations = $(stats.pumpIterations)")
        println("fwIterations = $(stats.fwIterations)")
        println("perturbCount = $(stats.perturbCount)")
        println("restartCount = $(stats.restartCount)")
    end

    println("exitReason = $exitMsg")
end