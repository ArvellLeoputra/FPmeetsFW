include("src/dependencies.jl")

include("src/utils/config.jl")
include("src/utils/stats.jl")
include("src/utils/display.jl")

include("src/scip/setup.jl")
include("src/scip/queries.jl")
include("src/scip/submit.jl")

include("src/fw/lmoMOI.jl")
include("src/fw/lmoLPI.jl")
include("src/fw/driver.jl")

include("src/heur/pump.jl")
include("src/heur/fpfw.jl")
include("src/heur/runner.jl")

fileName, config = loadConfig(ARGS)

if config.enablePlot
    include("src/utils/plot.jl")
end

runInstance(fileName, config)