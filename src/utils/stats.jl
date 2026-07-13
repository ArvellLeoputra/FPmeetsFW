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
        stats.exitReason = :scip_optimal
    elseif status == SCIP.SCIP_STATUS_INFEASIBLE
        stats.exitReason = :scip_infeasible
    elseif status == SCIP.SCIP_STATUS_UNBOUNDED
        stats.exitReason = :scip_unbounded
    elseif status == SCIP.SCIP_STATUS_TIMELIMIT
        stats.exitReason = :scip_time_limit
    elseif status == SCIP.SCIP_STATUS_NODELIMIT
        stats.exitReason = :scip_node_limit
    else
        stats.exitReason = :scip_unknown
    end

    stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)

    return stats
end

function printRunInfo(scip::SCIP.SCIPData)
    name = unsafe_string(SCIP.SCIPgetProbName(scip))
    varCount = SCIP.SCIPgetNOrigVars(scip)
    binCount = SCIP.SCIPgetNBinVars(scip)
    intCount = SCIP.SCIPgetNIntVars(scip)
    contCount = SCIP.SCIPgetNContVars(scip)

    printstyled("[run info]\n", color=:cyan)
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
    println("lineSearch = $(config.lineSearch)")
    println("timeLimit = $(config.timeLimit)")
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
    println("fwWarmStart = $(config.fwWarmStart ? "enabled" : "disabled")")
    println("lmoWarmStart = $(config.lmoWarmStart ? "enabled" : "disabled")")
    println("useSubMIP = $(config.useSubMIP ? "enabled" : "disabled")")
    println("presolve = $(config.presolve ? "enabled" : "disabled")")
    println("seed = $(config.seed)")
    println("enablePlot = $(config.enablePlot ? "enabled" : "disabled")")
    println("verbose = $(config.verbose ? "enabled" : "disabled")")
end

function printResults(stats::FPFWStats)
    exitMsg = if stats.exitReason == :time_limit
        "time limit reached"
    elseif stats.exitReason == :infeasible_fw
        "FW returned a point outside the feasible polytope (numerical error)"
    elseif stats.exitReason == :solution_found
        "integer feasible solution accepted by SCIP at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :rr_solution_found
        "integer feasible solution found by randomized rounding at iteration $(stats.pumpIterations)"
    elseif stats.exitReason == :scip_optimal
        "problem solved to optimality by SCIP before heuristic was called"
    elseif stats.exitReason == :scip_infeasible
        "LP is infeasible, thus MIP is infeasible"
    elseif stats.exitReason == :scip_time_limit
        "time limit reached before heuristic was called"
    elseif stats.exitReason == :scip_node_limit
        "SCIP node limit reached before heuristic was called"
    elseif stats.exitReason == :scip_unbounded
        "LP is unbounded, thus MIP is unbounded"
    else
        "unknown exit"
    end

    printstyled("[result]\n", color=:cyan)
    println("primalBound = $(round(stats.primalBound, digits=4))")
    println("dualBound = $(round(stats.dualBound, digits=4))")
    println("gap = $(@sprintf("%.2f %%", stats.gap * 100))")
    println("primalIntegral = $(round(stats.primalIntegral, digits=4))")
    println("solFound = $(stats.solutionFound)")
    println("totalTime = $(round(stats.totalTime, digits=2))s")

    if !startswith(string(stats.exitReason), "scip_")
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