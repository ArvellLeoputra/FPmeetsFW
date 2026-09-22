"""
FPFWConfig
A struct that holds the configuration parameters for the FPFW heuristic.
The fields can be set via a configuration file or command-line arguments.
"""
@kwdef struct FPFWConfig
    # --- General ---
    "name of the run"
    runName::String = "fpfw"
    "time budget in seconds for the FPFW heuristic only"
    timeLimit::Float64 = 300.0
    "SCIP's own presolve/root-LP deadline; 0 reuses timeLimit"
    scipTimeLimit::Float64 = 0.0
    "apply presolving to the original problem before running the heuristic"
    presolve::Bool = true
    "random seed for reproducibility"
    seed::Int = 42
    "verbosity level (0: summary only, 1: pump table, 2: per-iteration debug, 3: iteration diagnostics)"
    verbose::Int = 0

    # --- Feasibility Pump ---
    "objective weight for objective feasibility pump"
    alpha::Float64 = DEF_ALPHA
    "geometric decay factor applied to alpha each pump iteration"
    alphaFactor::Float64 = DEF_ALPHA_FACTOR
    "pump rounding step: use a randomized threshold (Bertacco et al., 2007)"
    randRound::Bool = false
    "probabilistically round the current LP point x as a cheap feasibility attempt;
    independent of the pump's own rounding"
    randFeasCheck::Bool = false
    "WalkSAT perturb: when there aren't enough fractional variables to flip,
    pull extras from the support of currently-violated LP rows"
    walksatPerturb::Bool = true
    "fallback each iteration: when xRound isn't accepted, fix integers to xRound and
    solve the LP to recover feasible continuous values"
    useDive::Bool = false

    # --- Frank-Wolfe ---
    "norm used for the projection step"
    norm::Symbol = :manhattan
    "variant of the Frank-Wolfe algorithm"
    fwVariant::Symbol = :vanilla
    "maximum number of iterations for the Frank-Wolfe algorithm"
    fwMaxIterations::Int = 1
    "step size strategy for the Frank-Wolfe algorithm"
    fwStepSize::Symbol = :unitary
    "constant step size gamma used only when fwStepSize == :fixed (curvature-free, so
    safe with a non-smooth norm like manhattan, unlike :adaptive/:secant/:backtracking)"
    fixedStepSize::Float64 = 0.5
    "warm start the active set for the Frank-Wolfe algorithm"
    fwWarmStart::Bool = false
    "round from a weighted-random active-set atom instead of xFrac (needs fwWarmStart, fwVariant != :vanilla)"
    fwActiveSetRound::Bool = false
    "use the raw SCIP-LPI LMO and seed its basis from SCIP's LP basis (vs the MOI-based LMO)"
    lmoWarmStart::Bool = true
    "stop a FW projection early after this many non-improving dual_gap steps in a row; 0 disables"
    fwStallMaxIters::Int = DEF_FW_STALL_MAX_ITERS
    "relative dual_gap improvement needed to reset fwStallMaxIters's counter"
    fwMinImprovement::Float64 = DEF_FW_MIN_IMPROVEMENT
    "consecutive stable pump iterations before a smooth-norm variable's LMO bounds get locked; 0 disables"
    fwLockThreshold::Int = DEF_FW_LOCK_THRESHOLD
    "max fraction of ncols newly locked per pump iteration"
    fwLockMaxFrac::Float64 = DEF_FW_LOCK_MAX_FRAC
    "absolute ceiling on new locks per iteration, regardless of ncols"
    fwLockCap::Int = DEF_FW_LOCK_CAP
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

    "time tracking in seconds"    
    totalTime::Float64 = 0.0
    rootTime::Float64 = 0.0
    heurTime::Float64 = 0.0
    rrTime::Float64 = 0.0
    fwTime::Float64 = 0.0
    diveTime::Float64 = 0.0
    lmoSetupTime::Float64 = 0.0

    "iteration statistics"
    pumpIterations::Int = 0
    fwIterations::Int = 0
    lmoSimplexIters::Int = 0

    "FW dual-gap tracking: sum/count of the final dual_gap reached by each FW-projection call, for a per-solve mean"
    sumFwGap::Float64 = 0.0
    fwGapSamples::Int = 0
    "times xRound was rounded from a sampled active-set atom instead of xFrac"
    activeSetRoundCount::Int = 0

    "lock-feature statistics"
    lockVarsPeak::Int = 0
    lockInfeasibleRetries::Int = 0
    lockRejected::Int = 0
    lockCheckSimplexIters::Int = 0
    
    "times randFeasCheck! / diveSolve were called"
    rfcCalls::Int = 0
    diveCalls::Int = 0
    "times submitSolution(xRound) actually ran (vs pumpIterations — the gap is the repeatRound skip)"
    roundSubmits::Int = 0

    "perturbation and restart statistics"
    perturbCount::Int = 0
    restartCount::Int = 0
    "times the WalkSAT fallback fired inside perturb"
    walksatCount::Int = 0
    "pump iteration at the first successful perturb, 0 if none occurred"
    firstPerturbation::Int = 0

    "solution status"
    solutionFound::Bool = false
    exitReason::FPFWExitReason = NONE
end

"""
LPILMO
A struct that implements a linear minimization oracle (LMO) for the 
Frank-Wolfe algorithm using SCIP's LPI interface.
"""
@kwdef struct LPILMO <: FrankWolfe.LinearMinimizationOracle
    "LPI pointer"
    lpi::Ptr{SCIP.SCIP_LPI}
    "number of columns and rows in the LMO"
    ncols::Int32
    nrows::Int32
    "original LP row count, LB, and UB before the manhattan aux block"
    origNrows::Int32
    origLb::Vector{Float64}
    origUb::Vector{Float64}
    "verbosity level"
    verbose::Int
    "deadline to stop the current LP solve"
    deadline::Base.RefValue{Float64}
    "cumulative dual-simplex iterations"
    simplexIters::Base.RefValue{Int}
    "most recent viable basis seen right after a real solve with no locks active"
    fallbackCstat::Base.RefValue{Vector{Cint}}
    fallbackRstat::Base.RefValue{Vector{Cint}}
end

"""
LMODeadlineExceeded
Thrown when an LP solve runs past its deadline. Caller treats it as a clean time-limit exit.
"""
struct LMODeadlineExceeded <: Exception end

"""
LMOPrimalInfeasible
Thrown when the LP solve finds no feasible point — usually locked variables clashing.
Caller unlocks everything and retries once.
"""
struct LMOPrimalInfeasible <: Exception end

"""
LMOSolverError
Thrown when the solver itself errors (SCIP_LPERROR), not a real infeasibility.
Caller gives up on this projection cleanly instead of crashing.
"""
struct LMOSolverError <: Exception end

"""
StageState
Per-stage mutable state of the FPFW pump loop. Reset at stage transition by `transitionToStage2!`,
except `stage1NoImpr`, which is a stage-1-only signal and is never reset on transition.
"""
@kwdef mutable struct StageState
    "active stage: 1 (binaries only) or 2 (all integers)"
    stage::Int
    "active integer index set"
    activeIntIdx::Vector{Int}
    "activeIntIdx as a Set, for O(1) lookup"
    activeIntIdxSet::Set{Int} = Set(activeIntIdx)
    "active general integer index set"
    activeGIntIdx::Vector{Int}
    "average number of variable flips per perturbation/restart"
    avgFlips::Int = min(DEF_AVG_FLIPS, length(activeIntIdx))
    "iteration count of current stage"
    stageIter::Int = 0
    "hash for last rounded point"
    prevHash::UInt = UInt(0)
    "lowest FW projection objective (distance from the projected point to its rounding target)
    since the last reset (reset on transition, perturb, and restart)"
    bestProjObj::Float64 = Inf
    "iterations with no improvement in bestProjObj, reset on transition/perturb/restart
    and if there is significant improvement (DEF_MIN_IMPROVEMENT)"
    stagnationCount::Int = 0
    "perturbs since the last restart or stage start"
    consecutivePerturbs::Int = 0
    "stage-1 only: iterations since the last significant (DEF_MIN_IMPROVEMENT) gain vs closestDist
    (the stage's all-time best), so it survives perturb/restart. Forces the stage 1 -> 2 handoff."
    stage1NoImpr::Int = 0
    "lowest FW projection objective seen this stage, reset only on transition"
    closestDist::Float64 = Inf
    "FW projection that achieved closestDist (LP-feasible); the stage-2 FW start at the transition"
    closestFrac::Vector{Float64}
    "the rounding that projection was aimed at; the stage-2 rounding target on the first iteration"
    closestRound::Vector{Float64}
    "set by transitionToStage2!, cleared after one iteration: the next rounding targets closestRound"
    transitioned::Bool = false
    "consecutive stable pump iterations per variable (xRound[i] unchanged)"
    roundStableCount::Vector{Int}
    "variables currently locked (lb=ub=xRound[i]) in the LMO"
    lockedIdx::Set{Int} = Set{Int}()
end

"""
LPInfo
A struct that holds the LP data for the FPFW heuristic,
received from getLPInfo function in src/scip/queries.jl
"""
@kwdef struct LPInfo
    "columns and rows of the LP, as SCIP pointers"
    lpCols::Vector{Ptr{SCIP.SCIP_COL}}
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}}
    "coefficient per LP column"
    objCoeffs::Vector{Float64}
    "objective function scale (euclidean norm of objCoeffs)"
    objScale::Float64
    "column index dictionary: SCIP column pointer -> index in lpCols"
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int}
    "indices of binary variables"
    binIdx::Vector{Int}
    "binIdx as a Set, for O(1) lookup"
    binSet::Set{Int}
    "indices of general integer variables"
    gIntIdx::Vector{Int}
    "indices of all integer variables"
    intIdx::Vector{Int}
    "number of columns and rows in the LP"
    ncols::Int32
    nrows::Int32
    "initial solution vector (LP-feasible)"
    initSol::Vector{Float64}
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
    "restart's general-integer flip budget: escalates across consecutive restarts, decays
    geometrically per iteration without one, capped at gIntFlipBudgetCap"
    gIntFlipBudget::Int = 0
    "pumpIterations value as of the last restart, used to compute the geometric decay above"
    lastRestart::Int = 0
    "cap on gIntFlipBudget, set once from the problem's general-integer count"
    gIntFlipBudgetCap::Int = 0
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