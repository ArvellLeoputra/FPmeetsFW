# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using Random
using Printf
import MathOptInterface
const MOI = MathOptInterface

struct FPFWConfig
    norm::Symbol
    fwVariant::Symbol
    lineSearch::Symbol
    randRound::Bool
    randFeasCheck::Bool
    warmStart::Bool
    presolve::Bool
    seed::Int
end

mutable struct FPFWStats
    primalBound::Float64
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
    exitReason::Symbol  # :none, :time_limit, :restart_limit, :infeasible_fw, :solution_found, :rr_solution_found, :scip_optimal, :scip_infeasible, :scip_time_limit, :scip_node_limit, :scip_unbounded, :scip_unknown

    FPFWStats() = new(Inf, Inf, Inf, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, false, :none)
end

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
    config::FPFWConfig
    startTime::Float64
    stats::FPFWStats
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
const DEF_INT_TOLERANCE = 1e-6
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

# Randomized rounding feasibility check parameters
const DEF_RAND_FEAS_ITER_LIMIT = 100

# Debug mode: set to true to print detailed step-by-step output
const DEBUG_VERBOSE = false

# Valid options for FPFWConfig fields
const VALID_NORMS = (:euclidean, :manhattan, :smoothManhattan)
const VALID_FW_VARIANTS = (:vanilla, :away, :blended_pairwise, :blended)
const VALID_LINE_SEARCHES = (:agnostic, :backtracking, :secant, :adaptive, :unitary)