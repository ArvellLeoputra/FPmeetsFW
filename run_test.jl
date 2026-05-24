include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")
include("fw_utils.jl")

function mps_test_model(fileName::String, config::FPFWConfig, globalStartTime::Float64)
    model = minimal_setup(presolve=DEF_PRESOLVE)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    SCIP.SCIPreadProb(scip, fileName, C_NULL)

    # Count original binary variables from MPS before any solving/presolving
    variableCount = SCIP.SCIPgetNOrigVars(scip)
    originalVars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetOrigVars(scip), variableCount)
    binaryCount = sum(SCIP.SCIPvarGetType(originalVars[j]) == SCIP.SCIP_VARTYPE_BINARY for j in 1:variableCount)
    integerCount = sum(SCIP.SCIPvarGetType(originalVars[j]) == SCIP.SCIP_VARTYPE_INTEGER for j in 1:variableCount)
    continuousCount = variableCount - binaryCount - integerCount

    printstyled("[run info]\n", color=:cyan)
    println("instance = $(splitext(basename(fileName))[1])")
    println("totalVars = $variableCount")
    println("binaryVars = $binaryCount")
    println("integerVars = $integerCount")
    println("continuousVars = $continuousCount")
    println("presolve = $(DEF_PRESOLVE ? "enabled" : "disabled")")

    printstyled("[FPFW configs]\n", color=:cyan)
    println("projectionNorm = $(config.projectionNorm)")
    println("fwVariant = $(config.fwVariant)")
    println("lineSearch = $(config.lineSearch)")
    println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
    println("warmStart = $(config.warmStart ? "enabled" : "disabled")")
    println("roundThreshold = $DEF_ROUNDING_THRESHOLD")

    heur = FPFWHeuristic(
        0,
        nothing,
        config,
        globalStartTime
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
        stats = FPFWStats()
        stats.totalTime = time() - globalStartTime
        stats.primalBound = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
        stats.gap = Float64(SCIP.SCIPgetGap(scip))
        stats.solutionFound = SCIP.SCIPgetNSols(scip) > 0                                                                                                                                                                         
        stats.exitReason = stats.solutionFound ? :scip_solved : :scip_time_limit                                                                                                                                                
        print_heuristic_summary(stats)
    end

    return nothing
end

if length(ARGS) < 1
    error("Usage: julia --project run_test.jl <fileName.mps> [euclidean|manhattan|smooth_manhattan] [vanilla|away|blended_pairwise|blended] [agnostic|backtracking|secant|adaptive]")
end

fileName = ARGS[1]

if !isfile(fileName)
    error("File not found: $fileName")
end

projectionNorm = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :euclidean

validNorms = (:euclidean, :manhattan, :smooth_manhattan)
if projectionNorm ∉ validNorms
    error("Invalid projection norm: $projectionNorm. Must be one of: $validNorms")
end

fwVariant = length(ARGS) >= 3 ? Symbol(ARGS[3]) : DEF_FW_VARIANT

validVariants = (:vanilla, :away, :blended_pairwise, :blended)
if fwVariant ∉ validVariants
    error("Invalid FW variant: $fwVariant. Must be one of: $validVariants")
end

lineSearch = length(ARGS) >= 4 ? Symbol(ARGS[4]) : DEF_LINE_SEARCH

validLineSearches = (:agnostic, :backtracking, :secant, :adaptive, :unitary)
if lineSearch ∉ validLineSearches
    error("Invalid line search: $lineSearch. Must be one of: $validLineSearches")
end

if projectionNorm == :manhattan && lineSearch ∈ (:adaptive, :secant)
    error("manhattan norm requires a smooth objective — use agnostic or backtracking line search instead")
end

config = FPFWConfig(projectionNorm, fwVariant, lineSearch, DEF_RAND_ROUND, DEF_WARM_START)
startTime = time()

if DEF_RANDOM_SEED !== nothing
    Random.seed!(DEF_RANDOM_SEED)
end

mps_test_model(fileName, config, startTime)
