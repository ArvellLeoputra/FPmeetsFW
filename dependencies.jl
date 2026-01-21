# Centralized Dependencies for FPFW Heuristic
module FPFWDependencies

using SCIP
using FrankWolfe

# MathOptInterface for building LMO
import MathOptInterface
const MOI = MathOptInterface

export MOI, FPFWHeuristic, FPFWStats
export DEF_TOLERANCE, DEF_FW_MAX_ITER, DEF_FP_MAX_ITER
export DEF_TIME_LIMIT, DEF_PERTURB_FRACTION, DEF_MAX_RESTARTS, DEF_RANDOM_SEED

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
    projection_norm::Symbol
end

mutable struct FPFWStats
    total_time::Float64
    fw_time::Float64
    fw_iterations::Int
    fp_iterations::Int
    solution_found::Bool
    final_objective::Union{Nothing, Float64}
    fw_calls::Int
    memory_used::Union{Nothing, Float64}

    FPFWStats() = new(0.0, 0.0, 0, 0, false, nothing, 0, 0.0)
end

# Default tolerance for feasibility/integrality checks
const DEF_TOLERANCE = 1e-6

# Iteration parameters
const DEF_FW_MAX_ITER = 100
const DEF_FP_MAX_ITER = 1000

# Time limit
const DEF_TIME_LIMIT = 3600.0

# Perturbation parameters
const DEF_PERTURB_FRACTION = 0.2      # Fraction of binary vars to flip when perturbing
const DEF_MAX_RESTARTS = 50           # Maximum number of restarts after cycles

# Random seed for reproducibility
const DEF_RANDOM_SEED = 42

end

using .FPFWDependencies

using JuMP
using SCIP
using FrankWolfe
using GLPK
using Random
using Printf
import MathOptInterface
