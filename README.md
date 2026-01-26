# FPmeetsFW

A primal heuristic for Mixed-Integer Programming (MIP) that combines the Feasibility Pump algorithm with Frank-Wolfe projection.

## Overview

The Feasibility Pump (FP) is a heuristic for finding feasible solutions to MIPs. It alternates between:
1. **Rounding** an LP-feasible solution to obtain an integer point
2. **Projecting** the rounded point back onto the LP relaxation

This implementation replaces the standard LP projection with **Frank-Wolfe** optimization, which projects onto the LP feasible region by minimizing a distance function.

## Algorithm

```
x <- LP solution
repeat:
    x_round <- round(x)           # Round binary variables
    x <- FrankWolfe(min ||x - x_round||, s.t. x is LP feasible)  # Project back
until x is integral or max iterations
```

## Features

- **Two projection norms**:
  - `euclidean`: Minimizes L2 distance (smooth, standard FP)
  - `manhattan`: Minimizes L1 distance (non-smooth)

- **Cycle detection**: Detects when the same rounded solution is visited twice

- **Perturbation**: Escapes cycles by randomly flipping binary variables in the rounding target

- **SCIP integration**: Runs as a SCIP primal heuristic

## Usage

```bash
julia --project run_test.jl <instance.mps> [euclidean|manhattan]
```

Examples:
```bash
julia --project run_test.jl ./testcase/test1.mps euclidean
julia --project run_test.jl ./testcase/test1.mps manhattan
```

## Parameters

Configurable in `dependencies.jl`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `DEF_FW_MAX_ITER` | 100 | Max Frank-Wolfe iterations per projection |
| `DEF_FP_MAX_ITER` | 1000 | Max Feasibility Pump iterations |
| `DEF_SCIP_TIME_LIMIT` | 3600s | Time limit |
| `DEF_FW_TIME_LIMIT` | 300s | Time limit |
| `DEF_PERTURB_FRACTION` | 0.2 | Fraction of binary vars to flip on perturbation |
| `DEF_MAX_RESTARTS` | 50 | Max restarts after cycle detection |
| `DEF_ROUNDING_THRESHOLD` | 0.47 | Rounding threshold |

## Dependencies

- [SCIP.jl](https://github.com/scipopt/SCIP.jl) - MIP solver
- [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) - Frank-Wolfe algorithm
- JuMP, GLPK

## File Structure

```
├── run_test.jl      # Entry point
├── dependencies.jl  # Parameters and type definitions
├── fpfwheur.jl      # Main FPFW heuristic implementation
├── lmo_builder.jl   # Builds Linear Minimization Oracle from SCIP LP
├── helper.jl        # Utility functions
├── scip_setup.jl    # SCIP configuration
└── testcase/        # Test instances (.mps files)
```
