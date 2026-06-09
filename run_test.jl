include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")
include("fw_utils.jl")

function mps_test_model(fileName::String, config::FPFWConfig, globalStartTime::Float64)
    Random.seed!(config.seed)
    model = minimal_setup(presolve=config.presolve)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner.scip[]

    GC.@preserve model begin
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
        println("presolve = $(config.presolve ? "enabled" : "disabled")")

        printstyled("[FPFW configs]\n", color=:cyan)
        println("norm = $(config.norm)")
        println("fwVariant = $(config.fwVariant)")
        println("lineSearch = $(config.lineSearch)")
        println("randomizedRounding = $(config.randRound ? "enabled" : "disabled")")
        println("randomizedFeasibilityCheck = $(config.randFeasCheck ? "enabled" : "disabled")")
        println("warmStart = $(config.warmStart ? "enabled" : "disabled")")
        println("seed = $(config.seed)")

        lpData = LPData(0, 0, Ptr{SCIP.SCIP_COL}[], Ptr{SCIP.SCIP_ROW}[], Dict{Ptr{SCIP.SCIP_COL}, Int}(), Int[], Int[], Int[], Float64[], nothing)
        solved = presolve_and_clone(backend, lpData)

        if solved
            stats = FPFWStats()
            stats.totalTime = time() - globalStartTime
            stats.primalBound = SCIP.SCIPgetNSols(scip) > 0 ? Float64(SCIP.SCIPgetPrimalbound(scip)) : nothing
            stats.gap = Float64(SCIP.SCIPgetGap(scip))
            stats.solutionFound = SCIP.SCIPgetNSols(scip) > 0
            stats.exitReason = :scip_presolved
            printHeurSummary(stats)
            return nothing
        end

        runFPFW(scip, lpData, config, globalStartTime)
    end
end

if length(ARGS) < 1
    error("Usage: julia --project run_test.jl <fileName.mps> [norm=] [fwVariant=] [fwLineSearch=] [randomizedRounding=] [randomFeasibilityCheck=] [warmStart=] [presolve=] [seed=]")
end

fileName = ARGS[1]

if !isfile(fileName)
    error("File not found: $fileName")
end

kwargs = Dict{String, String}()

# Load fpfw.cfg as defaults if it exists
if isfile("fpfw.cfg")
    for line in readlines("fpfw.cfg")
        line = strip(line)
        if !isempty(line) && !startswith(line, "#") && contains(line, "=")
            key, value = split(line, "=", limit=2)
            kwargs[key] = value
        end
    end
end

# Command-line args override cfg values
for arg in ARGS[2:end]
    if contains(arg, "=")
        key, value = split(arg, "=")
        kwargs[key] = value
    end
end

norm = Symbol(get(kwargs, "norm", string(DEF_NORM)))
validNorms = (:euclidean, :manhattan, :smoothManhattan)
if norm ∉ validNorms
    error("Invalid norm: $norm. Must be one of: $validNorms")
end

fwVariant = Symbol(get(kwargs, "fwVariant", string(DEF_FW_VARIANT)))
validVariants = (:vanilla, :away, :blended_pairwise, :blended)
if fwVariant ∉ validVariants
    error("Invalid fwVariant: $fwVariant. Must be one of: $validVariants")
end

lineSearch = Symbol(get(kwargs, "fwLineSearch", string(DEF_LINE_SEARCH)))
validLineSearches = (:agnostic, :backtracking, :secant, :adaptive, :unitary)
if lineSearch ∉ validLineSearches
    error("Invalid fwLineSearch: $lineSearch. Must be one of: $validLineSearches")
end

if norm == :manhattan && lineSearch ∈ (:adaptive, :secant)
    error("manhattan norm requires a smooth objective — use agnostic or backtracking instead")
end

randRound = parse(Bool, get(kwargs, "randomizedRounding", string(DEF_RAND_ROUND)))
randFeasCheck = parse(Bool, get(kwargs, "randomFeasibilityCheck", string(DEF_RAND_FEAS_CHECK)))
warmStart = parse(Bool, get(kwargs, "warmStart", string(DEF_WARM_START)))
presolve = parse(Bool, get(kwargs, "presolve", string(DEF_PRESOLVE)))
seed = parse(Int, get(kwargs, "seed", string(DEF_RANDOM_SEED)))

config = FPFWConfig(norm, fwVariant, lineSearch, randRound, randFeasCheck, warmStart, presolve, seed)
startTime = time()

mps_test_model(fileName, config, startTime)