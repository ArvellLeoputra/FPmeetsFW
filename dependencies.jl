# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using Random
using Printf
using Dates
import MathOptInterface
const MOI = MathOptInterface

struct FPFWConfig
    norm::Symbol
    fwVariant::Symbol
    fwMaxIterations::Int
    lineSearch::Symbol
    timeLimit::Float64
    randRound::Bool
    randFeasCheck::Bool
    fwWarmStart::Bool
    lmoWarmStart::Bool
    useSubMIP::Bool
    presolve::Bool
    seed::Int
    enablePlot::Bool
    verbose::Bool
end

mutable struct FPFWStats
    primalBound::Float64
    dualBound::Float64
    gap::Float64
    primalIntegral::Float64
    primalEvents::Vector{Tuple{Float64, Float64}}
    totalTime::Float64
    heurTime::Float64
    rrTime::Float64
    fwTime::Float64
    pumpIterations::Int
    fwIterations::Int
    perturbCount::Int
    restartCount::Int
    solutionFound::Bool
    exitReason::Symbol  # :none, :time_limit, :infeasible_fw, :solution_found, :rr_solution_found, :scip_optimal, :scip_infeasible, :scip_time_limit, :scip_node_limit, :scip_unbounded, :scip_unknown

    FPFWStats() = new(Inf, Inf, Inf, 0.0, Tuple{Float64,Float64}[], 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, false, :none)
end

struct BaseLMO <: FrankWolfe.LinearMinimizationOracle
    lpi::Ptr{SCIP.SCIP_LPI}
    ncols::Int32
    nrows::Int32
    verbose::Bool
end

mutable struct FPFWRunData
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO, BaseLMO}
    stats::FPFWStats

    FPFWRunData() = new(0, nothing, FPFWStats())
end

mutable struct FPFWHeuristic <: SCIP.Heuristic
    config::FPFWConfig
    data::FPFWRunData
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

# Perturbation and restart parameters
const DEF_MAX_STAGNATION = 3                # Maximum number of iterations without improvement before perturbing
const DEF_STAGNATION_RESTART_THRESHOLD = 3  # Maximum number of stagnations before restarting
const DEF_BIGM = 1e9                        # Big M constant for cycle-breaking perturbations
const DEF_BIGBIGM = 1e15                    # Bigbig M constant for perturbations

# Randomized rounding feasibility check parameters
const DEF_RAND_FEAS_ITER_LIMIT = 100

# Valid options for FPFWConfig fields
const VALID_NORMS = (:euclidean, :manhattan, :smoothManhattan)
const VALID_FW_VARIANTS = (:vanilla, :away, :blended_pairwise, :blended)
const VALID_LINE_SEARCHES = (:agnostic, :backtracking, :secant, :adaptive, :unitary)