function runInstance(fileName::String, config::FPFWConfig)
    startTime= time()
    Random.seed!(config.seed)

    model = minimal_setup(presolve=config.presolve)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, fileName, C_NULL)

    printRunInfo(scip, fileName)
    printConfigs(config)

    heur = FPFWHeuristic(0, nothing, config, startTime)
    SCIP.include_heuristic(
        backend,
        heur,
        name="FPFWHeuristic",
        priority=9999,
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )
    SCIP.SCIPsolve(scip)

    if heur.called == 0
        stats = FPFWStats()
        stats.totalTime = time() - startTime
        stats.solutionFound = SCIP.SCIPgetNSols(scip) > 0

        if stats.solutionFound
            stats.primalBound = Float64(SCIP.SCIPgetPrimalbound(scip))
            stats.gap = Float64(SCIP.SCIPgetGap(scip))
        end

        status = SCIP.SCIPgetStatus(scip)

        if status == SCIP.SCIP_STATUS_OPTIMAL || status == SCIP.SCIP_STATUS_INFEASIBLE
            stats.dualBound = Float64(SCIP.SCIPgetDualbound(scip))
        end

        if status == SCIP.SCIP_STATUS_OPTIMAL
            stats.exitReason = :scip_optimal  # presolve or root LP solved the problem
        elseif status == SCIP.SCIP_STATUS_INFEASIBLE
            stats.exitReason = :scip_infeasible  # infeasibility detected by presolve or root LP
        elseif status == SCIP.SCIP_STATUS_TIMELIMIT
            stats.exitReason = :scip_time_limit  # time limit hit during presolve or root LP
        elseif status == SCIP.SCIP_STATUS_NODELIMIT
            stats.exitReason = :scip_node_limit  # node limit hit during presolve or root LP
        elseif status == SCIP.SCIP_STATUS_UNBOUNDED
            stats.exitReason = :scip_unbounded  # unboundedness detected by presolve or root LP
        else
            stats.exitReason = :scip_unknown
        end

        printSCIPSummary(stats)
    end

    return nothing
end