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

function buildConfig(params::Dict{String, String})
    inputs = ["norm", "fwVariant", "fwLineSearch", "timeLimit", "randomizedRounding", "randomFeasibilityCheck", "fwWarmStart", "lmoWarmStart", "useSubMIP", "presolve", "seed"]
    for input in inputs
        if !haskey(params, input)
            error("Missing key '$input' in fpfw.cfg or CLI args")
        end
    end

    norm = Symbol(params["norm"])
    fwVariant = Symbol(params["fwVariant"])
    fwLineSearch = Symbol(params["fwLineSearch"])
    timeLimit = parse(Float64, params["timeLimit"])
    randRound = parse(Bool, params["randomizedRounding"])
    randFeasCheck = parse(Bool, params["randomFeasibilityCheck"])
    fwWarmStart = parse(Bool, params["fwWarmStart"])
    lmoWarmStart = parse(Bool, params["lmoWarmStart"])
    useSubMIP = parse(Bool, params["useSubMIP"])
    presolve = parse(Bool, params["presolve"])
    seed = parse(Int, params["seed"])

    if norm ∉ VALID_NORMS
        error("Invalid norm: $norm. Must be one of: $VALID_NORMS")
    end

    if fwVariant ∉ VALID_FW_VARIANTS
        error("Invalid fwVariant: $fwVariant. Must be one of: $VALID_FW_VARIANTS")
    end

    if fwLineSearch ∉ VALID_LINE_SEARCHES
        error("Invalid fwLineSearch: $fwLineSearch. Must be one of: $VALID_LINE_SEARCHES")
    end

    if norm == :manhattan && fwLineSearch ∈ (:adaptive, :secant, :backtracking)
        error("manhattan norm requires a smooth objective — use agnostic or unitary instead")
    end

    return FPFWConfig(norm, fwVariant, fwLineSearch, timeLimit, randRound, randFeasCheck, fwWarmStart, lmoWarmStart, useSubMIP, presolve, seed)
end

function loadConfig(args::Vector{String})
    if length(args) < 1
        error("Usage: julia main.jl <mps_file> [key=value ...]")
    end

    if !isfile(args[1])
        error("File not found: $(args[1])")
    end
    fileName = args[1]

    cfgPath = joinpath(@__DIR__, "fpfw.cfg")
    params = parseCfgFile(cfgPath)

    for arg in args[2:end]
        if contains(arg, "=")
            key, value = split(arg, "=", limit=2)
            params[strip(key)] = strip(value)
        end
    end

    return fileName, buildConfig(params)
end
