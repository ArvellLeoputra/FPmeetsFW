include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")

function mps_test_model(filename::String, projection_norm::Symbol, rounding_threshold::Float64, fw_variant::Symbol, line_search::Symbol)
    model = minimal_setup(verbosity=3)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, filename, C_NULL)

    heur = FPFWHeuristic(0, nothing, projection_norm, rounding_threshold, fw_variant, line_search)
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
    error("Usage: julia --project run_test.jl <filename.mps> [euclidean|manhattan] [threshold] [vanilla|away|blended_pairwise|blended] [agnostic|backtracking]")
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

rounding_threshold = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : DEF_ROUNDING_THRESHOLD

if rounding_threshold < 0.0 || rounding_threshold > 1.0
    error("Rounding threshold must be between 0.0 and 1.0, got: $rounding_threshold")
end

fw_variant = length(ARGS) >= 4 ? Symbol(ARGS[4]) : DEF_FW_VARIANT

valid_variants = (:vanilla, :away, :blended_pairwise, :blended)
if fw_variant ∉ valid_variants
    error("Invalid FW variant: $fw_variant. Must be one of: $valid_variants")
end

line_search = length(ARGS) >= 5 ? Symbol(ARGS[5]) : DEF_LINE_SEARCH

valid_line_searches = (:agnostic, :backtracking)
if line_search ∉ valid_line_searches
    error("Invalid line search: $line_search. Must be one of: $valid_line_searches")
end

println("\n" * "="^80)
println("Loading instance: $filename")
println("Projection norm: $projection_norm")
println("Rounding threshold: $rounding_threshold")
println("FW variant: $fw_variant")
println("Line search: $line_search")
println("="^80 * "\n")

start_time = time()
result = mps_test_model(filename, projection_norm, rounding_threshold, fw_variant, line_search)
total_time = time() - start_time

println("\n" * "="^80)
println("FINAL RESULT")
println("="^80)
println("Status:       $(result.status)")
println("Solutions:    $(result.nsols)")
println("Objective:    $(result.objective)")
println("Total time:   $(round(total_time, digits=2)) seconds")
println("="^80)