include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")

function mps_test_model(filename::String, projection_norm::Symbol, rounding_threshold::Float64, fw_variant::Symbol, line_search::Symbol, presolve::Bool, global_start_time::Float64)
    model = minimal_setup(presolve=presolve)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, filename, C_NULL)

    # Count original binary variables from MPS before any solving/presolving
    nvars_orig = SCIP.SCIPgetNOrigVars(scip)
    orig_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetOrigVars(scip), nvars_orig)
    n_orig_binary = sum(SCIP.SCIPvarGetType(orig_vars[j]) == SCIP.SCIP_VARTYPE_BINARY for j in 1:nvars_orig)

    println("="^80)
    println("RUN INFO")
    println("="^80)
    println("Instance:          $(basename(filename))")
    println("Total variables:   $nvars_orig")
    println("Binary variables:  $n_orig_binary")
    println("Projection norm:   $projection_norm")
    println("Rounding thresh:   $rounding_threshold")
    println("FW variant:        $fw_variant")
    println("Line search:       $line_search")
    println("Presolve:          $presolve")
    println("="^80)

    heur = FPFWHeuristic(n_orig_binary, 0, nothing, projection_norm, rounding_threshold, fw_variant, line_search, global_start_time)
    SCIP.include_heuristic(
        backend,
        heur,
        name="FPFWHeuristic",
        priority=9999,
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )

    SCIP.SCIPsolve(scip)

    if heur.called == 0
        total_time = time() - global_start_time
        objective = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
        gap = Float64(SCIP.SCIPgetGap(scip))
        stats = FPFWStats()
        stats.exit_reason = :scip_time_limit
        print_heuristic_summary(stats, total_time, objective, gap)
    end

    return nothing
end

if length(ARGS) < 1
    error("Usage: julia --project run_test.jl <filename.mps> [euclidean|manhattan|abssmooth] [threshold] [vanilla|away|blended_pairwise|blended] [agnostic|backtracking|secant|adaptive] [true|false]")
end

filename = ARGS[1]

if !isfile(filename)
    error("File not found: $filename")
end

projection_norm = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :euclidean

valid_norms = (:euclidean, :manhattan, :abssmooth)
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

valid_line_searches = (:agnostic, :backtracking, :secant, :adaptive)
if line_search ∉ valid_line_searches
    error("Invalid line search: $line_search. Must be one of: $valid_line_searches")
end

if projection_norm == :manhattan && line_search ∈ (:adaptive, :secant)
    error("manhattan norm requires a smooth objective — use agnostic or backtracking line search instead")
end

presolve = length(ARGS) >= 6 ? parse(Bool, ARGS[6]) : DEF_PRESOLVE

start_time = time()
mps_test_model(filename, projection_norm, rounding_threshold, fw_variant, line_search, presolve, start_time)
