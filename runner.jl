function runInstance(fileName::String, config::FPFWConfig)
    startTime= time()
    Random.seed!(config.seed)

    model = minimal_setup(presolve=config.presolve)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, fileName, C_NULL)

    printRunInfo(scip, fileName)
    printConfigs(config)

    heur = FPFWHeuristic(0, nothing, config, startTime, FPFWStats())
    SCIP.include_heuristic(
        backend,
        heur,
        name="FPFWHeuristic",
        priority=9999,
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )
    SCIP.SCIPsolve(scip)

    stats = buildStats(scip, startTime, heur)
    printResults(stats)

    return nothing
end