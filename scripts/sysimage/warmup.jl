# Precompile execution file for the FP-FW sysimage (see build.jl).
#
# Runs the pump end-to-end for every config in settings/, at both verbosity
# levels, over a small set of instances chosen so that between them every
# expensive code path is executed and therefore compiled into the image:
#
#   neos5      - solved at iteration 1 by FW projection (feasFWProj accept path)
#   gen-ip002  - solved at iteration 1 by direct rounding (feasRound accept path)
#   30n20b8    - hundreds of perturbations and restarts, stage-1 -> stage-2
#                transition (LMO rebuild), diving (with useDive), and the
#                iteration-/time-limit exit paths
#
# A held-out instance is used afterwards with --trace-compile to confirm that no
# method of the pump or its numerical dependencies is compiled during a real run.
#
# Env overrides:
#   FPFW_WARMUP_INSTANCES   comma-separated instance basenames (default: the set above)
#   FPFW_WARMUP_TIMELIMIT   per-run heuristic time cap in seconds (default: 30)

using FPmeetsFW

const SYSIMG_DIR   = @__DIR__
const FPFW_DIR     = normpath(joinpath(SYSIMG_DIR, "..", ".."))
const SETTINGS_DIR = joinpath(FPFW_DIR, "settings")
const INSTANCE_DIR = normpath(joinpath(FPFW_DIR, "..", "instances", "miplib_selected"))
const TLIM         = get(ENV, "FPFW_WARMUP_TIMELIMIT", "30")

const DEFAULT_INSTANCES = ["neos5", "gen-ip002", "30n20b8"]
const INSTANCES = split(get(ENV, "FPFW_WARMUP_INSTANCES", join(DEFAULT_INSTANCES, ",")), ",")

function resolve_instance(name)
    for ext in (".mps.gz", ".mps")
        p = joinpath(INSTANCE_DIR, name * ext)
        isfile(p) && return p
    end
    error("warmup instance not found: $name in $INSTANCE_DIR")
end

cfgs = sort(filter(f -> endswith(f, ".cfg"), readdir(SETTINGS_DIR)))
@info "sysimage warmup" instances = INSTANCES configs = cfgs timelimit = TLIM

for name in INSTANCES
    mps = resolve_instance(String(strip(name)))
    for cfg in cfgs, verbose in ("1", "2")
        cfgpath = joinpath(SETTINGS_DIR, cfg)
        args = String[mps, cfgpath, "timeLimit=$TLIM", "verbose=$verbose"]
        try
            fileName, config, _ = loadConfig(args)
            redirect_stdout(devnull) do
                runInstance(fileName, config)
            end
            @info "warmup ok" instance = name cfg verbose
        catch e
            @warn "warmup run errored (continuing)" instance = name cfg verbose exception = (e, catch_backtrace())
        end
    end
end

@info "sysimage warmup done"
