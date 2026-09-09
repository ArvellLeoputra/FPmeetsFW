# Central place for all tunable constants of the FPFW heuristic.

# Tolerances for feasibility/integrality checks and FW convergence
const DEF_INT_TOLERANCE = 1e-6
const DEF_FW_TOLERANCE = 1e-7

# Below this gradient magnitude the iterate is essentially at its rounding target, so a sign
# change is floating-point noise, not a real crossing (used by the debug grad-flip detector)
const DEF_GRAD_FLIP_TOL = 1e-3

# SCIP's own limits/time, set effectively unbounded so it never cuts the solve off early
const DEF_SCIP_TIME_LIMIT = 1000000

# Perturbation and restart parameters
const DEF_MAX_STAGNATION = 3       # Iterations without a real (>= threshold) projObj improvement before perturbing
const DEF_MIN_IMPROVEMENT = 0.10   # Relative improvement in projObj required to reset the stagnation counter
const DEF_MAX_PERTURBS = 10        # Cumulative perturbs in current stage before escalating to a restart
const DEF_AVG_FLIPS = 20           # Fixed perturbation flip count
const DEF_BIGM = 1e9               # Big M constant for cycle-breaking perturbations
const DEF_BIGBIGM = 1e15           # Bigbig M constant for perturbations

# Objective Feasibility Pump parameters
const DEF_ALPHA = 0.0  # initial objective weight
const DEF_ALPHA_FACTOR = 0.9  # geometric decay per iteration
const DEF_ALPHA_MIN = 1e-4  # alpha at or below this snaps to 0

# Staging parameters
const DEF_STAGE1_MAX_ITER = 10000
const DEF_STAGE2_MAX_ITER = 2000
const DEF_STAGE1_NOIMPR_LIMIT = DEF_STAGE1_MAX_ITER ÷ 10  # Stage 1 iteration limit without improvement in bestProjObj

# Global cap on total pump iterations
const DEF_MAX_PUMP_ITER = 12000

# Randomized rounding feasibility check parameters
const DEF_RAND_FEAS_ITER_LIMIT = 100

# Pump display formatting parameters
const DEF_MAX_INT_DIGITS = 7  # switch a float column to scientific notation beyond this many integer digits