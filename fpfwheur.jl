using SCIP

mutable struct FPFWHeuristic <: SCIP.Heuristic
    max_iter::Int
end

function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::SCIP.Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR}
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}

    @info("FWFP heuristic called")
    result = SCIP.SCIP_DIDNOTFIND
    return (SCIP.SCIP_OKAY, result)
end

optimizer = SCIP.Optimizer()

heuristic_storage = Dict{Any, Ptr{SCIP.SCIP_HEUR}}()

SCIP.include_heuristic(
    optimizer,
    FPFWHeuristic(100),
    name="FPFWHeuristic",
    description="Frank–Wolfe Feasibility Pump heuristic",
    dispchar='F',
    timing_mask = SCIP.SCIP_HEURTIMING_DURINGLPLOOP
)

scip = optimizer.inner
SCIP.@SCIP_CALL(SCIP.SCIPreadProb(scip, "gen-ip002.mps", C_NULL))
SCIP.@SCIP_CALL(SCIP.SCIPsolve(scip))