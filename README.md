# FPmeetsFW

A primal heuristic for Mixed-Integer Programming (MIP) that combines the Feasibility Pump algorithm with Frank-Wolfe projection.

## Overview

The Feasibility Pump (FP) is a heuristic for finding feasible solutions to MIPs. It alternates between:
1. **Rounding** an LP-feasible solution to obtain an integer point
2. **Projecting** the rounded point back onto the LP relaxation

This implementation replaces the standard LP projection with **Frank-Wolfe** optimization, which projects onto the LP feasible region by minimizing a distance function.

The **Frank-Wolfe** (FW) algorithm (also known as the conditional gradient method) is an iterative first-order method for constrained optimization. Instead of projecting the gradient step back onto the feasible set (as in projected gradient methods), each FW iteration solves a simpler **Linear Minimization Oracle (LMO)** — finding the vertex of the feasible region in the direction of the negative gradient. The next iterate is then a convex combination of the current point and this vertex, which keeps all iterates feasible by construction. This makes FW particularly well-suited for LP-feasible regions, where the LMO is just an LP solve.

## Algorithm

```
x <- LP solution
repeat:
    xRound <- round(x)           # Round binary variables
    x <- FrankWolfe(min ||x - xRound||, s.t. x is LP feasible)  # Project back
until x is integral or time limit / max iterations
```

## Features

- **Three projection norms**:
  - `euclidean`: Minimizes L2 distance (smooth, standard FP)
  - `manhattan`: Minimizes L1 distance (non-smooth)
  - `smoothManhattan`: Minimizes L1 distance with a smooth approximation

- **Four FW variants**: `vanilla`, `away`, `blended_pairwise`, `blended`

- **Five line search strategies**: `unitary`, `agnostic`, `backtracking`, `secant`, `adaptive`

- **Cycle detection**:
  - *Rounding cycle*: same rounded solution visited again → perturb rounding target

- **SCIP integration**: Runs as a SCIP primal heuristic at the root node

## Usage

```bash
julia --project main.jl <fileName.mps> [key=value ...]
```

All arguments are optional and fall back to `settings/fpfw.cfg`. The cfg file is required.

Examples:
```bash
julia --project main.jl ./testcase/test1.mps
julia --project main.jl ./testcase/test1.mps norm=euclidean seed=123
```

## Configuration

Settings are loaded in priority order: **command-line args > `settings/fpfw.cfg`**

`settings/fpfw.cfg` (run configuration):

| Key | Default | Description |
|-----|---------|-------------|
| `norm` | `manhattan` | Projection norm (`euclidean`, `manhattan`, `smoothManhattan`) |
| `fwVariant` | `vanilla` | FW variant (`vanilla`, `away`, `blended_pairwise`, `blended`) |
| `fwMaxIterations` | `1` | Max Frank-Wolfe iterations per FP projection step |
| `fwLineSearch` | `unitary` | Line search (`unitary`, `agnostic`, `backtracking`, `secant`, `adaptive`) |
| `timeLimit` | `300.0` | Total heuristic time limit (seconds) |
| `randomizedRounding` | `false` | Randomized rounding threshold |
| `randomFeasibilityCheck` | `false` | Randomized feasibility check |
| `fwWarmStart` | `true` | Warm-start FW active set across FP iterations |
| `lmoWarmStart` | `true` | Warm-start the LMO's LP basis across FP iterations |
| `useSubMIP` | `false` | Fall back to fixing integer vars and resolving the LP when rounding fails |
| `presolve` | `true` | Enable SCIP presolving |
| `seed` | `42` | Random seed for reproducibility |
| `enablePlot` | `false` | Plot FW iterates (2D instances only) |
| `verbose` | `false` | Print detailed per-iteration output |

`src/dependencies.jl` (algorithm constants):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DEF_INT_TOLERANCE` | `1e-6` | Tolerance for feasibility/integrality checks |
| `DEF_FW_TOLERANCE` | `1e-7` | FW convergence tolerance (duality gap) |
| `DEF_MAX_STAGNATION` | `3` | Max iterations without improvement before perturbing |
| `DEF_STAGNATION_RESTART_THRESHOLD` | `3` | Max stagnation-triggered perturbs before restarting |
| `DEF_BIGM` / `DEF_BIGBIGM` | `1e9` / `1e15` | Big-M constants used when perturbing general integer variables |
| `DEF_RAND_FEAS_ITER_LIMIT` | `100` | Max attempts for randomized feasibility rounding |
| `DEF_MAX_INT_DIGITS` | `7` | Switch a pump-table float column to scientific notation beyond this many integer digits |

## Dependencies

- [SCIP.jl](https://github.com/scipopt/SCIP.jl) - MIP solver
- [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) - Frank-Wolfe algorithm
- JuMP

## File Structure

```
├── main.jl                  # Entry point — parses args and calls runInstance
├── settings/
│   └── fpfw.cfg             # Run configuration
├── src/
│   ├── dependencies.jl      # Structs, constants, and valid option sets
│   ├── fw/
│   │   ├── driver.jl        # Builds FW objective/gradient/line-search, dispatches to FrankWolfe.jl variants
│   │   ├── lmoLPI.jl        # LMO backed by SCIP's raw LPI (default, warm-startable)
│   │   └── lmoMOI.jl        # LMO backed by MathOptInterface (cold-start fallback)
│   ├── heur/
│   │   ├── fpfwheur.jl      # Main FPFW heuristic implementation (SCIP.find_primal_solution)
│   │   ├── pump.jl          # Rounding, cycle-detection hashing, perturb/restart mechanics
│   │   └── runner.jl        # runInstance — sets up and runs the solver
│   ├── scip/
│   │   ├── setup.jl         # JuMP-level SCIP solver configuration
│   │   ├── queries.jl       # Read-only SCIP/LP state queries
│   │   └── submit.jl        # Solution submission and sub-MIP solving (mutates SCIP state)
│   └── utils/
│       ├── config.jl        # Configuration loading and validation
│       ├── stats.jl         # Stats aggregation and console reporting
│       ├── display.jl       # Pump-table formatting
│       └── plot.jl          # 2D iterate plotting (optional)
└── testcase/                # Sample .mps instances
```