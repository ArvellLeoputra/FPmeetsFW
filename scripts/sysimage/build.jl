# Build the FP-FW sysimage: bakes SCIP / FrankWolfe / JuMP / MOI + FPmeetsFW and
# every specialization the warmup exercises into one image, so sweep tasks start
# with zero JIT and near-zero package load time.
#
# Usage:
#   julia scripts/sysimage/build.jl [output_path]
#
# Default output: <FPmeetsFW>/fpfw_sysimage.so  (or $FPFW_SYSIMAGE)
#
# Run the pump with it:
#   julia --sysimage <output_path> --project=<FPmeetsFW> main.jl <instance> <cfg> ...
#   submit_sweep.sh -S <output_path> ...    # wires --sysimage into the array job
#
# Rebuild whenever src/, settings/, Project.toml or a dependency changes.

import Pkg

# PackageCompiler is a build-time tool - keep it out of FPmeetsFW's own deps.
try
    @eval import PackageCompiler
catch
    @info "installing PackageCompiler into the default environment"
    Pkg.activate(; temp = false)          # default (@vX.Y) environment
    Pkg.add("PackageCompiler")
    @eval import PackageCompiler
end
using PackageCompiler

const SYSIMG_DIR = @__DIR__
const FPFW_DIR   = normpath(joinpath(SYSIMG_DIR, "..", ".."))
const OUT = abspath(
    !isempty(ARGS)               ? ARGS[1] :
    haskey(ENV, "FPFW_SYSIMAGE") ? ENV["FPFW_SYSIMAGE"] :
    joinpath(FPFW_DIR, "fpfw_sysimage.so")
)

@info "building FP-FW sysimage" project = FPFW_DIR output = OUT

PackageCompiler.create_sysimage(
    ["FPmeetsFW"];
    project                   = FPFW_DIR,
    sysimage_path             = OUT,
    precompile_execution_file = joinpath(SYSIMG_DIR, "warmup.jl"),
    cpu_target                = "generic",   # portable across cluster nodes;
                                             # set "native" only if building on a
                                             # node identical to the compute nodes
)

@info "sysimage written" OUT
println(OUT)
