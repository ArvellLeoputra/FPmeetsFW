using SCIP
using JuMP
using FrankWolfe
using GLPK
import MathOptInterface
const MOI = MathOptInterface

include("scipsetup.jl")

include("fpfwheur.jl")

function mps_test_model()
    model = minimal_setup()
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner
    
    SCIP.SCIPreadProb(scip, "filename.mps", C_NULL)  # Insert MPS file here

    heur = FPFWHeuristic(0, nothing)
    SCIP.include_heuristic(
        backend, 
        heur,
        name="FPFWHeuristic", 
        priority=9999, 
        timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
    )

    SCIP.SCIPsolve(scip)
end

mps_test_model()