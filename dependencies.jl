# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using Random
using Printf
import MathOptInterface
const MOI = MathOptInterface

struct FPFWConfig
    projectionNorm::Symbol
    fwVariant::Symbol
    lineSearch::Symbol
    randRound::Bool
    warmStart::Bool
end

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
    config::FPFWConfig 
    globalStartTime::Float64
end

mutable struct FPFWStats
    primalBound::Union{Float64, Nothing}
    dualBound::Float64
    gap::Float64
    totalTime::Float64
    heurTime::Float64
    rrTime::Float64
    fwTime::Float64
    pumpIterations::Int
    fwIterations::Int
    restartCount::Int
    solutionFound::Bool
    exitReason::Symbol  # :none, :time_limit, :restart_limit, :infeasible_fw, :solution_found, :rr_solution_found, :solution_rejected, :scip_time_limit, :scip_solved

    FPFWStats() = new(nothing, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, false, :none)
end

mutable struct PumpDisplayColumn
    name::String
    width::Int
    decimals::Int
end

mutable struct PumpDisplay
    column::Vector{PumpDisplayColumn}
end

# Default tolerance for feasibility/integrality checks and FW convergence
const DEF_TOLERANCE = 1e-6
const DEF_FW_TOLERANCE = 1e-7

# Iteration parameters
const DEF_FW_MAX_ITER = 1

# Time limit
const DEF_GLOBAL_TIME_LIMIT = 480.0
const DEF_SCIP_TIME_LIMIT = 300.0

# FW escape check: check if FW escapes its rounding point
const DEF_FW_ESCAPE = false

# Perturbation parameters
const DEF_PERTURB_FRACTION = 0.2   # Fraction of binary vars to flip when cycle detected
const DEF_MAX_RESTARTS = 1000      # Maximum number of restarts before giving up
const DEF_MAX_STAGNATION = 3       # Maximum number of iterations without improvement before perturbing
const DEF_BIGM = 1e9               # Big M constant for cycle-breaking perturbations
const DEF_BIGBIGM = 1e15           # Bigbig M constant for perturbations

# Random seed for reproducibility; set to nothing to disable
const DEF_RANDOM_SEED::Union{Nothing, Int} = 42

# Rounding threshold for deciding when to round fractional solutions
const DEF_ROUNDING_THRESHOLD = 0.5

# Randomized rounding parameters
const DEF_RAND_ROUND = true
const DEF_RR_TIME_LIMIT = 3.0

# Warm-starting away/blended variants with the previous iteration's active set
const DEF_WARM_START = true

# Frank-Wolfe variant: :vanilla, :away, :blended_pairwise, :blended
const DEF_FW_VARIANT = :vanilla

# Frank-Wolfe line search: :unitary, :agnostic, :backtracking, :secant, :adaptive
const DEF_LINE_SEARCH = :unitary

# Determine whether presolve on or off
const DEF_PRESOLVE = true

# Debug mode: set to true to print detailed step-by-step output
const DEBUG_VERBOSE = false