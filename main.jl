include("dependencies.jl")
include("config.jl")
include("scip_setup.jl")
include("helper.jl")
include("lmo_builder.jl")
include("fw_utils.jl")
include("fpfwheur.jl")
include("runner.jl")

fileName, config = loadConfig(ARGS)

if config.enablePlot
    include("plot_utils.jl")
end

runInstance(fileName, config)