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

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
end

# Default tolerance for feasibility/integrality checks
const DEF_TOLERANCE = 1e-6

# Default Frank-Wolfe parameters
const DEF_FW_MAX_ITER = 10
const DEF_FP_MAX_ITER = 100