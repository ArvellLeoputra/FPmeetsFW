# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using Random
using Printf
import MathOptInterface
const MOI = MathOptInterface

mutable struct FPFWHeuristic <: SCIP.Heuristic
    binary_vars::Int
    general_integer_vars::Int
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
    projection_norm::Symbol
    rounding_threshold::Float64
    fw_variant::Symbol      # :vanilla, :away, :blended_pairwise, :blended
    line_search::Symbol     # :agnostic, :backtracking, :secant, :adaptive
    global_start_time::Float64
end

mutable struct FPFWStats
    heur_time::Float64
    rr_time::Float64
    fw_time::Float64
    fw_iterations::Int
    fp_iterations::Int
    restarts::Int
    solution_found::Bool
    exit_reason::Symbol  # :none, :time_limit, :restart_limit, :infeasible_fw, :solution_found, :rr_solution_found, :solution_rejected, :scip_time_limit, :scip_solved
    iter_found_solution::Union{Nothing, Int}
    lp_objective::Float64
    final_objective::Union{Nothing, Float64}

    FPFWStats() = new(0.0, 0.0, 0.0, 0, 0, 0, false, :none, nothing, 0.0, nothing)
end

# Default tolerance for feasibility/integrality checks and FW convergence
const DEF_TOLERANCE = 1e-6
const DEF_FW_TOLERANCE = 1e-7

# Iteration parameters
const DEF_FW_MAX_ITER = 1000

# Time limit
const DEF_GLOBAL_TIME_LIMIT = 480.0
const DEF_SCIP_TIME_LIMIT = 300.0

# Perturbation parameters
const DEF_PERTURB_FRACTION = 0.2   # Fraction of binary vars to flip when cycle detected
const DEF_MAX_RESTARTS = 1000  # Maximum number of restarts before giving up
const DEF_BIGM = 1e9  # Big M constant for cycle-breaking perturbations
const DEF_BIGBIGM = 1e15  # Bigbig M constant for perturbations

# Random seed for reproducibility; set to nothing to disable
const DEF_RANDOM_SEED::Union{Nothing, Int} = 42

# Rounding threshold for deciding when to round fractional solutions;
const DEF_ROUNDING_THRESHOLD = 0.5

# Randomized rounding pre-check: n_attempts = div(n_integers, DEF_RAND_ROUND_DIVISOR)
const DEF_RAND_ROUND_DIVISOR = 2

# Frank-Wolfe variant: :vanilla, :away, :blended_pairwise, :blended
const DEF_FW_VARIANT = :vanilla

# Frank-Wolfe line search: :agnostic, :backtracking, :secant, :adaptive
const DEF_LINE_SEARCH = :agnostic

# Debug mode: set to true to print detailed step-by-step output
const DEBUG_VERBOSE = false

# Determine whether presolve on or off
const DEF_PRESOLVE = false