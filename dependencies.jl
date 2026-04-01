# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using GLPK
using Random
using Printf
import MathOptInterface
const MOI = MathOptInterface

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
    projection_norm::Symbol
    rounding_threshold::Float64
    fw_variant::Symbol      # :vanilla, :away, :blended_pairwise, :blended
    line_search::Symbol     # :agnostic, :backtracking
end

mutable struct FPFWStats
    total_time::Float64
    fw_time::Float64
    fw_iterations::Int
    fp_iterations::Int
    solution_found::Bool
    final_objective::Union{Nothing, Float64}

    FPFWStats() = new(0.0, 0.0, 0, 0, false, nothing)
end

# Default tolerance for feasibility/integrality checks
const DEF_TOLERANCE = 1e-6

# Iteration parameters
const DEF_FW_MAX_ITER = 100
const DEF_FP_MAX_ITER = 1000

# Time limit
const DEF_SCIP_TIME_LIMIT = 600.0
const DEF_FW_TIME_LIMIT = 300.0

# Perturbation parameters
const DEF_PERTURB_FRACTION = 0.2      # Fraction of binary vars to flip when perturbing
const DEF_MAX_RESTARTS = 50           # Maximum number of restarts after cycles

# Random seed for reproducibility; set to nothing to disable
const DEF_RANDOM_SEED::Union{Nothing, Int} = 42

# Rounding threshold for deciding when to round fractional solutions;
# For basic rounding, set to 0.5
const DEF_ROUNDING_THRESHOLD = 0.47

# Frank-Wolfe variant: :vanilla, :away, :blended_pairwise, :blended
const DEF_FW_VARIANT = :vanilla

# Frank-Wolfe line search: :agnostic, :backtracking
const DEF_LINE_SEARCH = :agnostic

# Debug mode: set to true to print detailed step-by-step output
const DEBUG_VERBOSE = false
