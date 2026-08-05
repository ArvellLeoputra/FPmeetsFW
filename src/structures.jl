"""
FPFWConfig
A struct that holds the configuration parameters for the FPFW heuristic.
The fields can be set via a configuration file or command-line arguments.
"""
@kwdef struct FPFWConfig
    "name of the run"
    runName::String = "fpfw"
    "norm used for the projection step"
    norm::Symbol = :manhattan
    "variant of the Frank-Wolfe algorithm"
    fwVariant::Symbol = :vanilla
    "maximum number of iterations for the Frank-Wolfe algorithm"
    fwMaxIterations::Int = 1
    "step size strategy for the Frank-Wolfe algorithm"
    fwStepSize::Symbol = :unitary
    "total solving time budget in seconds (SCIP + heuristic)"
    timeLimit::Float64 = 300.0
    "pump rounding step: use a randomized threshold (Bertacco et al., 2007)"
    randRound::Bool = false
    "probabilistically round the current LP point x as a cheap feasibility attempt;
    independent of the pump's own rounding"
    randFeasCheck::Bool = false
    "warm start the active set for the Frank-Wolfe algorithm"
    fwWarmStart::Bool = false
    "use the raw SCIP-LPI LMO and seed its basis from SCIP's LP basis (vs the MOI-based LMO)"
    lmoWarmStart::Bool = true
    "fallback each iteration: when xRound isn't accepted, fix integers to xRound and
    solve the LP to recover feasible continuous values"
    useDive::Bool = false
    "apply presolving to the original problem before running the heuristic"
    presolve::Bool = true
    "random seed for reproducibility"
    seed::Int = 42
    "verbosity level (0: summary only, 1: pump table, 2: per-iteration debug, 3: iteration diagnostics)"
    verbose::Int = 0
end

"""
FPFWExitReason
An enum for possible termination status of the FPFW heuristic
"""
@enum FPFWExitReason begin
    NONE
    # Heuristic-specific exit reasons
    TIME_LIMIT
    ITER_LIMIT
    INFEASIBLE_FW
    SOLUTION_ROUND     # xRound accepted directly
    SOLUTION_DIVE      # integers fixed to xRound, LP solved for continuous
    SOLUTION_FWPROJ    # FW-projected point accepted
    SOLUTION_RR        # randomized feasibility check
    # SCIP_* mean SCIP finished before the heuristic ran
    SCIP_OPTIMAL
    SCIP_INFEASIBLE
    SCIP_TIME_LIMIT
    SCIP_NODE_LIMIT
    SCIP_UNBOUNDED
    SCIP_UNKNOWN
end

"""
FPFWStats
A struct that holds the statistics of a run of the FPFW heuristic.
Used for logging and reporting purposes.
"""
@kwdef mutable struct FPFWStats
    "objective bounds and final gap"
    primalBound::Float64 = Inf
    dualBound::Float64 = Inf
    gap::Float64 = Inf

    "primal integral tracking"
    primalIntegral::Float64 = 0.0
    primalEvents::Vector{Tuple{Float64, Float64}} = Tuple{Float64, Float64}[]

    "time breakdown in seconds: total solve, heuristic wall time, randomized feasibility-check time, FW-solve time"
    totalTime::Float64 = 0.0
    heurTime::Float64 = 0.0
    rrTime::Float64 = 0.0
    fwTime::Float64 = 0.0

    "iteration statistics"
    pumpIterations::Int = 0
    fwIterations::Int = 0

    "perturbation and restart statistics"
    perturbCount::Int = 0
    restartCount::Int = 0

    "solution status"
    solutionFound::Bool = false
    exitReason::FPFWExitReason = NONE
end

"""
LPILMO
A struct that implements a linear minimization oracle (LMO) for the 
Frank-Wolfe algorithm using SCIP's LPI interface.
"""
struct LPILMO <: FrankWolfe.LinearMinimizationOracle
    lpi::Ptr{SCIP.SCIP_LPI}
    ncols::Int32
    nrows::Int32      # total rows (original LP rows + manhattan aux rows, if any)
    origNrows::Int32  # original LP row count, before the manhattan aux block
    verbose::Int
end

"""
FPFWRunData
A struct that holds the runtime data for a run of the FPFW heuristic.
"""
@kwdef mutable struct FPFWRunData
    "number of times SCIP has invoked the heuristic"
    called::Int64 = 0
    "the Frank-Wolfe linear minimization oracle: MathOptLMO (MOI-based) or LPILMO (raw SCIP-LPI)"
    lmo::Union{Nothing, FrankWolfe.MathOptLMO, LPILMO} = nothing
    "run statistics, filled during the run and read out afterwards for reporting"
    stats::FPFWStats = FPFWStats()
    "MOI-path only: per general-integer, the two |x - xRound| constraint handles, rewritten each
    iteration with the current rounding target. The LPI path needs no equivalent — LPIupdateRounding!
    addresses the aux rows by their integer position (appended after the original nrows rows) instead."
    auxConstraintRefs::Vector{Tuple{MOI.ConstraintIndex, MOI.ConstraintIndex}} = Tuple{MOI.ConstraintIndex, MOI.ConstraintIndex}[]
end

"""
FPFWHeuristic
A struct that implements the FPFW heuristic as a SCIP plugin.
"""
mutable struct FPFWHeuristic <: SCIP.Heuristic
    config::FPFWConfig
    data::FPFWRunData
end