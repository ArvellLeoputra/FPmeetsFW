include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")

function mps_test_model(filename::String, projection_norm::Symbol)
    model = minimal_setup(verbosity=3)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, filename, C_NULL)

    heur = FPFWHeuristic(0, nothing, projection_norm)
    SCIP.include_heuristic(
        backend,
        heur,
        name="FPFWHeuristic",
        priority=9999,
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )

    SCIP.SCIPsolve(scip)

    # Return solution info
    status = SCIP.SCIPgetStatus(scip)
    nsols = SCIP.SCIPgetNSols(scip)
    obj = nsols > 0 ? SCIP.SCIPgetPrimalbound(scip) : nothing

    return (status=status, nsols=nsols, objective=obj)
end

if length(ARGS) < 1
    error("Usage: julia --project run_test.jl <filename.mps> [euclidean|manhattan]")
end

filename = ARGS[1]

if !isfile(filename)
    error("File not found: $filename")
end

projection_norm = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :euclidean

valid_norms = (:euclidean, :manhattan)
if projection_norm ∉ valid_norms
    error("Invalid projection norm: $projection_norm. Must be one of: $valid_norms")
end

println("\n" * "="^80)
println("Loading instance: $filename")
println("Projection norm: $projection_norm")
println("="^80 * "\n")

start_time = time()
result = mps_test_model(filename, projection_norm)
total_time = time() - start_time

println("\n" * "="^80)
println("FINAL RESULT")
println("="^80)
println("Status:       $(result.status)")
println("Solutions:    $(result.nsols)")
println("Objective:    $(result.objective)")
println("Total time:   $(round(total_time, digits=2)) seconds")
println("="^80)