include("src/FPmeetsFW.jl")
using .FPmeetsFW

fileName, config, resultsDir = loadConfig(ARGS)
runInstance(fileName, config, resultsDir)