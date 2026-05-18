include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")
include("fw_utils.jl")

function mps_test_model(filename::String, config::FPFWConfig, global_start_time::Float64)
    model = minimal_setup(presolve=DEF_PRESOLVE)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, filename, C_NULL)

    # Count original binary variables from MPS before any solving/presolving
    nvars_orig = SCIP.SCIPgetNOrigVars(scip)
    orig_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetOrigVars(scip), nvars_orig)
    n_orig_binary = sum(SCIP.SCIPvarGetType(orig_vars[j]) == SCIP.SCIP_VARTYPE_BINARY for j in 1:nvars_orig)
    n_orig_integer = sum(SCIP.SCIPvarGetType(orig_vars[j]) == SCIP.SCIP_VARTYPE_INTEGER for j in 1:nvars_orig)
    n_orig_continuous = nvars_orig - n_orig_binary - n_orig_integer

    println("="^80)
    println("RUN INFO")
    println("="^80)
    println("Instance:               $(basename(filename))")
    println("Total variables:        $nvars_orig")
    println("Binary variables:       $n_orig_binary")
    println("G integer variables:    $n_orig_integer")
    println("Continuous variables:   $n_orig_continuous")
    println("Presolve:               $(DEF_PRESOLVE ? "enabled" : "disabled")")
    println("="^80)

    println()
    println("="^80)
    println("HEURISTIC CONFIGURATION")
    println("="^80)
    println("Projection norm:        $(config.projection_norm)")
    println("FW variant:             $(config.fw_variant)")
    println("Line search:            $(config.line_search)")
    println("Randomized rounding:    $(config.rand_round ? "enabled" : "disabled")")
    println("Warm start:             $(config.warm_start ? "enabled" : "disabled")")
    println("Rounding threshold:     $DEF_ROUNDING_THRESHOLD")
    println("="^80)

    heur = FPFWHeuristic(
        0,
        nothing,
        config,
        global_start_time
    )

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

        if SCIP.SCIPgetNSols(scip) > 0
            stats.solution_found = true
            stats.exit_reason = :scip_solved
        else
            stats.exit_reason = :scip_time_limit
        end
        
        print_heuristic_summary(stats, total_time, objective, gap)
    end

    return nothing
end

if length(ARGS) < 1
    error("Usage: julia --project run_test.jl <filename.mps> [euclidean|manhattan|smooth_manhattan] [vanilla|away|blended_pairwise|blended] [agnostic|backtracking|secant|adaptive] [rand_round=true|false] [warm_start=true|false]")
end

filename = ARGS[1]

if !isfile(filename)
    error("File not found: $filename")
end

projection_norm = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :euclidean

valid_norms = (:euclidean, :manhattan, :smooth_manhattan)
if projection_norm ∉ valid_norms
    error("Invalid projection norm: $projection_norm. Must be one of: $valid_norms")
end

fw_variant = length(ARGS) >= 3 ? Symbol(ARGS[3]) : DEF_FW_VARIANT

valid_variants = (:vanilla, :away, :blended_pairwise, :blended)
if fw_variant ∉ valid_variants
    error("Invalid FW variant: $fw_variant. Must be one of: $valid_variants")
end

line_search = length(ARGS) >= 4 ? Symbol(ARGS[4]) : DEF_LINE_SEARCH

valid_line_searches = (:agnostic, :backtracking, :secant, :adaptive, :unitary)
if line_search ∉ valid_line_searches
    error("Invalid line search: $line_search. Must be one of: $valid_line_searches")
end

if projection_norm == :manhattan && line_search ∈ (:adaptive, :secant)
    error("manhattan norm requires a smooth objective — use agnostic or backtracking line search instead")
end

rand_round = length(ARGS) >= 5 ? parse(Bool, ARGS[5]) : DEF_RAND_ROUND
warm_start  = length(ARGS) >= 6 ? parse(Bool, ARGS[6]) : DEF_WARM_START

if warm_start && fw_variant == :vanilla
    @warn "warm_start is enabled but fw_variant=:vanilla does not support warm starting — warm start will be ignored"
end

config = FPFWConfig(projection_norm, fw_variant, line_search, rand_round, warm_start)
start_time = time()

if DEF_RANDOM_SEED !== nothing
    Random.seed!(DEF_RANDOM_SEED)
end

mps_test_model(filename, config, start_time)
