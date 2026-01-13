include("dependencies.jl")
include("scip_setup.jl")
include("helper.jl")
include("fpfwheur.jl")
include("lmo_builder.jl")

function mps_test_model(filename::String)
    model = minimal_setup()
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner
    
    SCIP.SCIPreadProb(scip, filename, C_NULL)  # Insert MPS file here

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

if length(ARGS) < 1
    error("No .mps file supplied. Usage:\n    julia script.jl <filename.mps>")
end

filename = ARGS[1]
println("Loading instance: $filename")

mps_test_model(filename)