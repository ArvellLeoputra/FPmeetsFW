function runInstance(fileName::String, config::FPFWConfig)
    printstyled("Timestamp: ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"), "\n", color=:cyan)

    globalStartTime = time()
    Random.seed!(config.seed)

    model = minimal_setup(time_limit=config.timeLimit, presolve=config.presolve, verbosity=config.verbose)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, fileName, C_NULL)

    printRunInfo(scip, config)
    printConfigs(config)

    heur = FPFWHeuristic(config, FPFWRunData())
    SCIP.include_heuristic(
        backend,
        heur,
        name="FPFWHeuristic",
        priority=9999,
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )
    SCIP.SCIPsolve(scip)

    totalTime = timeElapsed(globalStartTime)

    stats = buildStats(scip, heur, totalTime)
    printResults(stats)

    return nothing
end