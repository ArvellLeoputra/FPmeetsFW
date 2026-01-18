# Centralized Dependencies for FPFW Heuristic
using JuMP
using SCIP
using FrankWolfe
using GLPK

# MathOptInterface for building LMO
import MathOptInterface
import MathOptInterface:
    add_variables,
    add_constraint,
    copy_to,
    ScalarAffineFunction,
    ScalarAffineTerm,
    GreaterThan,
    LessThan

const MOI = MathOptInterface

using Printf

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
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

# Default Frank-Wolfe parameters
const DEF_FW_MAX_ITER = 100
const DEF_FP_MAX_ITER = 1000

# Time limit in seconds (5 minutes)
const DEF_TIME_LIMIT = 300.0
