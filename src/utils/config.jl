# Valid options for FPFWConfig fields
const VALID_NORMS = (:euclidean, :manhattan, :smoothManhattan)
const VALID_FW_VARIANTS = (:vanilla, :away, :blended_pairwise, :blended)
const VALID_FW_STEP_SIZES = (:agnostic, :backtracking, :secant, :adaptive, :unitary)

function parseCfgFile(path::String)
    if !isfile(path)
        error("Configuration file 'fpfw.cfg' not found: $path")
    end

    params = Dict{String, String}()
    for line in readlines(path)
        line = strip(line)
        if !isempty(line) && !startswith(line, "#") && contains(line, "=")
            key, value = split(line, "=", limit=2)
            params[strip(key)] = strip(value)
        end
    end
    return params
end

function buildFPFWConfig(params::Dict{String, String})
    inputs = [
        "runName",
        "norm",
        "fwVariant",
        "fwMaxIterations",
        "fwStepSize",
        "timeLimit",
        "randomizedRounding",
        "randomFeasibilityCheck",
        "fwWarmStart",
        "lmoWarmStart",
        "useDive",
        "presolve",
        "seed",
        "verbose"
    ]
    
    for input in inputs
        if !haskey(params, input)
            error("Missing key '$input' in fpfw.cfg or CLI args")
        end
    end

    runName = params["runName"]
    norm = Symbol(params["norm"])
    fwVariant = Symbol(params["fwVariant"])
    fwMaxIterations = parse(Int, params["fwMaxIterations"])
    fwStepSize = Symbol(params["fwStepSize"])
    timeLimit = parse(Float64, params["timeLimit"])
    randRound = parse(Bool, params["randomizedRounding"])
    randFeasCheck = parse(Bool, params["randomFeasibilityCheck"])
    fwWarmStart = parse(Bool, params["fwWarmStart"])
    lmoWarmStart = parse(Bool, params["lmoWarmStart"])
    useDive = parse(Bool, params["useDive"])
    presolve = parse(Bool, params["presolve"])
    seed = parse(Int, params["seed"])
    verbose = parse(Int, params["verbose"])

    if verbose < 0 || verbose > 3
        error("Invalid verbose: $verbose. Must be 0, 1, 2, or 3")
    end

    if norm ∉ VALID_NORMS
        error("Invalid norm: $norm. Must be one of: $VALID_NORMS")
    end

    if fwVariant ∉ VALID_FW_VARIANTS
        error("Invalid fwVariant: $fwVariant. Must be one of: $VALID_FW_VARIANTS")
    end

    if fwStepSize ∉ VALID_FW_STEP_SIZES
        error("Invalid fwStepSize: $fwStepSize. Must be one of: $VALID_FW_STEP_SIZES")
    end

    # Avoid manhattan norm with step sizes that require a smooth objective
    if norm == :manhattan && fwStepSize ∈ (:adaptive, :secant, :backtracking)
        error("manhattan norm requires a smooth objective — use agnostic or unitary instead")
    end

    # Avoid warm-starting the active set for manhattan norm
    if norm == :manhattan && fwWarmStart
        error("manhattan norm does not support fwWarmStart — set fwWarmStart=false")
    end

    if fwMaxIterations < 0
        error("Invalid fwMaxIterations: $fwMaxIterations. Must be a positive integer")
    end

    return FPFWConfig(
        runName = runName,
        norm = norm,
        fwVariant = fwVariant,
        fwMaxIterations = fwMaxIterations,
        fwStepSize = fwStepSize,
        timeLimit = timeLimit,
        randRound = randRound,
        randFeasCheck = randFeasCheck,
        fwWarmStart = fwWarmStart,
        lmoWarmStart = lmoWarmStart,
        useDive = useDive,
        presolve = presolve,
        seed = seed,
        verbose = verbose,
    )
end

function loadConfig(args::Vector{String})
    if length(args) < 1
        error("Usage: julia main.jl <mps_file> [key=value ...]")
    end

    # Load MPS file
    if !isfile(args[1])
        error("File not found: $(args[1])")
    end
    fileName = args[1]

    # Load configuration file
    settingsDir = joinpath(@__DIR__, "..", "..", "settings")
    cliArgsIdx = 2
    
    if length(args) >= 2 && !contains(args[2], "=")
        cfgArg = args[2]
        if isfile(cfgArg)  # absolute path
            cfgPath = cfgArg
        elseif isfile(joinpath(settingsDir, cfgArg))  # relative path
            cfgPath = joinpath(settingsDir, cfgArg)
        else
            error("Configuration file not found: $cfgArg")
        end
        cliArgsIdx = 3
    else
        cfgPath = joinpath(settingsDir, "fpfw.cfg")
    end

    params = parseCfgFile(cfgPath)

    # Override parameters with CLI args
    for arg in args[cliArgsIdx:end]
        if contains(arg, "=")
            key, value = split(arg, "=", limit=2)
            params[strip(key)] = strip(value)
        end
    end

    # resultsDir is an optional CLI-only key: the directory to write results.json into.
    resultsDir = get(params, "resultsDir", "")

    return fileName, buildFPFWConfig(params), resultsDir
end