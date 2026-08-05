module FPmeetsFW

using JuMP
using SCIP
using FrankWolfe
using Random
using Printf
using Dates
import MathOptInterface
const MOI = MathOptInterface

include("structures.jl")
include("constants.jl")

include("utils/config.jl")
include("utils/stats.jl")
include("utils/display.jl")

include("scip/setup.jl")
include("scip/queries.jl")
include("scip/submit.jl")

include("fw/lmoMOI.jl")
include("fw/lmoLPI.jl")
include("fw/driver.jl")

include("heur/pump.jl")
include("heur/fpfw.jl")
include("heur/runner.jl")

export loadConfig, runInstance

end