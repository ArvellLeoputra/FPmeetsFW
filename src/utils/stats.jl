function timeElapsed(start::Float64)
    return time() - start
end

# Instance name (no directory, no .mps/.mps.gz suffix) from an instance file path.
function instanceName(fileName::String)
    base = basename(fileName)
    for suf in (".mps.gz", ".mps")
        if endswith(base, suf)
            return base[1:end-length(suf)]
        end
    end
    return base
end

# instance name -> known optimal/best-known objective, from SOLU_FILE's "=opt="/"=best=" lines.
function loadKnownOptimal(fileName::String)
    name = instanceName(fileName)
    isfile(SOLU_FILE) || return NaN
    for line in eachline(SOLU_FILE)
        parts = split(line)
        length(parts) >= 3 || continue
        status = strip(parts[1], '=')
        (status == "opt" || status == "best") || continue
        if parts[2] == name
            return parse(Float64, parts[3])
        end
    end
    return NaN
end

# Primal gap vs a known optimal/best-known reference:
#   0 if both are 0; |x*-x~|/|x*| if same-sign and nonzero; Inf otherwise (incl. not found).
function computeGap(knownOptimal::Float64, primalBound::Float64, solutionFound::Bool)
    isnan(knownOptimal) && return NaN  # no reference available for this instance
    !solutionFound && return Inf
    if knownOptimal == 0.0 && primalBound == 0.0
        return 0.0
    elseif knownOptimal * primalBound > 0
        return abs(knownOptimal - primalBound) / abs(knownOptimal)
    else
        return Inf
    end
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

function buildStats(scip::SCIP.SCIPData, heur::FPFWHeuristic, totalTime::Float64, fileName::String)
    stats = heur.data.stats
    if heur.data.called > 0
        stats.totalTime = totalTime
        stats.primalIntegral = computePrimalIntegral(stats.primalEvents, totalTime)
        stats.knownOptimal = loadKnownOptimal(fileName)
        stats.gap = computeGap(stats.knownOptimal, stats.primalBound, stats.solutionFound)
        return stats
    end

    stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))

    stats.totalTime = totalTime
    stats.rootTime = totalTime
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
    stats.knownOptimal = loadKnownOptimal(fileName)
    stats.gap = computeGap(stats.knownOptimal, stats.primalBound, stats.solutionFound)

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
    println("readTime = $(round(SCIP.SCIPgetReadingTime(scip); digits=4))s")
end

function printConfigs(config::FPFWConfig)
    printstyled("[FPFW configs]\n", color=:cyan)
    println("norm = $(config.norm)")
    println("fwVariant = $(config.fwVariant)")
    println("fwMaxIterations = $(config.fwMaxIterations)")
    println("fwStepSize = $(config.fwStepSize)")
    if config.fwStepSize == :fixed
        println("fixedStepSize = $(config.fixedStepSize)")
    end
    println("fwStallMaxIters = $(config.fwStallMaxIters)")
    if config.fwStallMaxIters > 0
        println("fwMinImprovement = $(config.fwMinImprovement)")
    end
    println("fwLockThreshold = $(config.fwLockThreshold)")
    if config.fwLockThreshold > 0
        println("fwLockMaxFrac = $(config.fwLockMaxFrac)")
        println("fwLockCap = $(config.fwLockCap)")
    end
    println("timeLimit = $(config.timeLimit)")
    if config.scipTimeLimit > 0
        println("scipTimeLimit = $(config.scipTimeLimit)")
    end
    println("alpha = $(config.alpha)")
    println("alphaFactor = $(config.alphaFactor)")
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("walksatPerturb = $(config.walksatPerturb ? "enabled" : "disabled")")
    println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
    println("fwWarmStart = $(config.fwWarmStart ? "enabled" : "disabled")")
    println("fwActiveSetRound = $(config.fwActiveSetRound ? "enabled" : "disabled")")
    println("fwAsrShakes = $(config.fwAsrShakes)")
    println("lmoWarmStart = $(config.lmoWarmStart ? "enabled" : "disabled")")
    println("useDive = $(config.useDive ? "enabled" : "disabled")")
    println("exactCycleCheck = $(config.exactCycleCheck ? "enabled" : "disabled")")
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
    println("knownOptimal = $(isnan(stats.knownOptimal) ? "n/a" : round(stats.knownOptimal, digits=4))")
    gapStr = isnan(stats.gap) ? "n/a (no reference)" : isinf(stats.gap) ? "inf" : @sprintf("%.2f %%", stats.gap * 100)
    println("gap = $gapStr")
    println("primalIntegral = $(round(stats.primalIntegral, digits=4))")
    println("solFound = $(stats.solutionFound)")
    println("totalTime = $(round(stats.totalTime, digits=4))s")
    println("rootTime = $(round(stats.rootTime, digits=4))s")
    println("exitReason = $exitMsg")

    if !startswith(string(stats.exitReason), "SCIP_")
        printstyled("[pump stats]\n", color=:cyan)
        println("heurTime = $(round(stats.heurTime, digits=4))s")
        println("lmoSetupTime = $(round(stats.lmoSetupTime, digits=4))s")
        println("fwTime = $(round(stats.fwTime, digits=4))s")
        println("randRoundTime = $(round(stats.rrTime, digits=4))s")
        println("diveTime = $(round(stats.diveTime, digits=4))s")
        println("pumpIterations = $(stats.pumpIterations)")
        println("rfcCalls = $(stats.rfcCalls)")
        println("diveCalls = $(stats.diveCalls)")
        println("roundSubmits = $(stats.roundSubmits)")
        println("activeSetRoundCount = $(stats.activeSetRoundCount)")
        println("asrShakeCount = $(stats.asrShakeCount)")
        println("fwIterations = $(stats.fwIterations)")
        avgFwGapStr = stats.fwGapSamples > 0 ? @sprintf("%.6g", stats.sumFwGap / stats.fwGapSamples) : "n/a"
        println("avgFwGap = $avgFwGapStr")
        println("lmoSimplexIters = $(stats.lmoSimplexIters)")
        println("lockVarsPeak = $(stats.lockVarsPeak)")
        println("lockInfeasibleRetries = $(stats.lockInfeasibleRetries)")
        println("lockRejected = $(stats.lockRejected)")
        println("lockCheckSimplexIters = $(stats.lockCheckSimplexIters)")
        println("perturbCount = $(stats.perturbCount)")
        println("restartCount = $(stats.restartCount)")
        println("walksatCount = $(stats.walksatCount)")
        println("firstPerturbation = $(stats.firstPerturbation)")
    end
end

function writeResults(stats::FPFWStats, config::FPFWConfig, fileName::String, resultsDir::String)
    mkpath(resultsDir)

    result = Dict(
        "instance" => basename(fileName),
        "runName" => config.runName,
        "primalBound" => stats.primalBound,
        "knownOptimal" => isnan(stats.knownOptimal) ? nothing : stats.knownOptimal,
        # gap vs knownOptimal (Definition 1.2, Mexi 2021), not SCIP's own dual bound.
        # "Inf" = no solution found (or the formula's own sign-mismatch case);
        # "NaN" = no known-optimal reference available for this instance.
        "gap" => isnan(stats.gap) ? "NaN" : (isinf(stats.gap) ? "Inf" : stats.gap),
        "primalIntegral" => stats.primalIntegral,
        "solutionFound" => stats.solutionFound,
        "totalTime" => stats.totalTime,
        "rootTime" => stats.rootTime,
        "heurTime" => stats.heurTime,
        "fwTime" => stats.fwTime,
        "randRoundTime" => stats.rrTime,
        "diveTime" => stats.diveTime,
        "lmoSetupTime" => stats.lmoSetupTime,
        "rfcCalls" => stats.rfcCalls,
        "diveCalls" => stats.diveCalls,
        "roundSubmits" => stats.roundSubmits,
        "activeSetRoundCount" => stats.activeSetRoundCount,
        "asrShakeCount" => stats.asrShakeCount,
        "pumpIterations" => stats.pumpIterations,
        "fwIterations" => stats.fwIterations,
        "avgFwGap" => stats.fwGapSamples > 0 ? stats.sumFwGap / stats.fwGapSamples : nothing,
        "lmoSimplexIters" => stats.lmoSimplexIters,
        "lockVarsPeak" => stats.lockVarsPeak,
        "lockInfeasibleRetries" => stats.lockInfeasibleRetries,
        "lockRejected" => stats.lockRejected,
        "lockCheckSimplexIters" => stats.lockCheckSimplexIters,
        "perturbCount" => stats.perturbCount,
        "restartCount" => stats.restartCount,
        "walksatCount" => stats.walksatCount,
        "firstPerturbation" => stats.firstPerturbation,
        "exitReason" => string(stats.exitReason),
    )

    open(joinpath(resultsDir, "results.json"), "w") do io
        JSON.print(io, result, 2)
    end
end