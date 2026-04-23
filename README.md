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
    x_round <- round(x)           # Round binary variables
    x <- FrankWolfe(min ||x - x_round||, s.t. x is LP feasible)  # Project back
until x is integral or time limit / max iterations
```

## Features

- **Three projection norms**:
  - `euclidean`: Minimizes L2 distance (smooth, standard FP)
  - `manhattan`: Minimizes L1 distance (non-smooth)
  - `abssmooth`: Minimizes L1 distance with a smooth approximation

- **Four FW variants**: `vanilla`, `away`, `blended_pairwise`, `blended`

- **Four line search strategies**: `agnostic`, `backtracking`, `secant`, `adaptive`

- **Cycle detection**:
  - *Rounding cycle*: same rounded solution visited again → perturb rounding target

- **SCIP integration**: Runs as a SCIP primal heuristic at the root node

## Usage

```bash
julia --project run_test.jl <instance.mps> [euclidean|manhattan|abssmooth] [threshold] [vanilla|away|blended_pairwise|blended] [agnostic|backtracking|secant|adaptive] [true|false]
```

Examples:
```bash
julia --project run_test.jl ./testcase/test1.mps euclidean 0.5 vanilla secant false
julia --project run_test.jl ./testcase/test1.mps manhattan 0.47 away agnostic true
```

## Parameters

Configurable in `dependencies.jl`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DEF_FW_MAX_ITER` | 1000 | Max Frank-Wolfe iterations per FP projection step |
| `DEF_GLOBAL_TIME_LIMIT` | 480s | Total wall-clock time limit for the heuristic |
| `DEF_SCIP_TIME_LIMIT` | 300s | SCIP solver time limit (for root LP solve) |
| `DEF_TOLERANCE` | 1e-6 | Tolerance for feasibility/integrality checks |
| `DEF_FW_TOLERANCE` | 1e-7 | FW convergence tolerance (duality gap) |
| `DEF_PERTURB_FRACTION` | 0.2 | Fraction of binary vars to flip on cycle restart |
| `DEF_MAX_RESTARTS` | 1000 | Max restarts after cycle detection |
| `DEF_ROUNDING_THRESHOLD` | 0.5 | Threshold for rounding fractional values to 1 |
| `DEF_FW_VARIANT` | `:vanilla` | FW variant (`:vanilla`, `:away`, `:blended_pairwise`, `:blended`) |
| `DEF_LINE_SEARCH` | `:agnostic` | Line search (`:agnostic`, `:backtracking`, `:secant`, `:adaptive`) |
| `DEF_RANDOM_SEED` | `42` | Random seed for reproducibility (`nothing` to disable) |
| `DEF_PRESOLVE` | `false` | Whether to enable SCIP presolving before the heuristic |
| `DEBUG_VERBOSE` | `false` | Print detailed per-iteration output |

## Dependencies

- [SCIP.jl](https://github.com/scipopt/SCIP.jl) - MIP solver
- [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) - Frank-Wolfe algorithm
- JuMP

## File Structure

```
├── run_test.jl      # Entry point
├── dependencies.jl  # Parameters and type definitions
├── fpfwheur.jl      # Main FPFW heuristic implementation
├── lmo_builder.jl   # Builds Linear Minimization Oracle from SCIP LP
├── helper.jl        # Utility functions
└── scip_setup.jl    # SCIP configuration
```
