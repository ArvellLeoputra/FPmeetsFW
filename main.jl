include("src/FPmeetsFW.jl")
using .FPmeetsFW

fileName, config = loadConfig(ARGS)
runInstance(fileName, config)