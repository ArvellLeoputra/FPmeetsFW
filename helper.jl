function print_heuristic_summary(
    stats::FPFWStats,
    total_time::Float64,
    objective::Union{Float64, Nothing},
    gap::Float64
)::Nothing

    exit_msg = if stats.exit_reason == :time_limit
        "global time limit $(DEF_GLOBAL_TIME_LIMIT)s reached"
    elseif stats.exit_reason == :restart_limit
        "FP cycled $(DEF_MAX_RESTARTS) times without progress"
    elseif stats.exit_reason == :infeasible_fw
        "FW returned a point outside the feasible polytope (numerical error)"
    elseif stats.exit_reason == :solution_found
        "integer feasible solution accepted by SCIP at iteration $(stats.fp_iterations)"
    elseif stats.exit_reason == :rr_solution_found
        "integer feasible solution found by randomized rounding at iteration $(stats.fp_iterations)"
    elseif stats.exit_reason == :solution_rejected
        "integer feasible solution found but rejected by SCIP"
    elseif stats.exit_reason == :scip_time_limit
        "SCIP time limit $(DEF_SCIP_TIME_LIMIT)s exceeded, heuristic never called"
    elseif stats.exit_reason == :scip_solved                                                                                                                                                                                          
        "problem solved by SCIP presolve/LP before heuristic was called"
    else
        "unknown exit"
    end

    heur_time_str = stats.heur_time == 0.0 ? "N/A" : "$(round(stats.heur_time, digits=2))s"  # 0.0 means heuristic never ran
    obj_str       = objective === nothing ? "N/A" : "$(round(objective, digits=4))"
    gap_str       = isinf(gap) || gap > 1e15 ? "Infinite" : @sprintf("%.2f %%", gap * 100)

    println("\n" * "="^80)
    println("FPFW HEURISTIC SUMMARY")
    println("="^80)
    println("Objective:         $obj_str")
    println("Gap:               $gap_str")
    println("Total time:        $(round(total_time, digits=2))s")
    println("Total heur time:   $heur_time_str")
    println("FW time:           $(round(stats.fw_time, digits=2))s")
    println("RR time:           $(round(stats.rr_time, digits=2))s")
    println("FP iterations:     $(stats.fp_iterations)")
    println("FW iterations:     $(stats.fw_iterations)")
    println("Restarts:          $(stats.restarts)")
    println("Solution found:    $(stats.solution_found)")
    println("Exit reason:       $exit_msg")
    println("="^80 * "\n")
end

function is_equal_values(
    x1::Float64,
    x2::Float64,
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    return abs(x1 - x2) <= tolerance
end

function is_lower_than(
    x1::Float64,
    x2::Float64,
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    return x1 < x2 - tolerance
end

function are_equal_vectors(
    integer_indices::Vector{Int},
    x1::Vector{Float64},
    x2::Vector{Float64},
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    for i in integer_indices
        if abs(x1[i] - x2[i]) > tolerance
            return false
        end
    end
    return true
end

function extract_lp_data(
    scip::Ptr{SCIP.SCIP_},
    nvars::Int32,
)
    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)
    col_to_idx = Dict(lp_cols[k] => k for k in 1:nvars)

    binary = Int[]
    integer = Int[]
    current_solution = zeros(SCIP.SCIP_Real, nvars)

    for j in 1:nvars
        var = SCIP.SCIPcolGetVar(lp_cols[j])
        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binary, j)
        elseif SCIP.SCIPvarIsIntegral(var) == SCIP.TRUE
            push!(integer, j)
        end
        current_solution[j] = SCIP.SCIPcolGetPrimsol(lp_cols[j])
    end

    return lp_cols, col_to_idx, binary, integer, current_solution
end

# Rounding function using custom threshold
function round_solution!(
    x_round::Vector{Float64},
    x::Vector{Float64},
    integer_indices::Vector{Int},
    threshold::Float64
)::Nothing

    for i in integer_indices
        x_round[i] = x[i] - floor(x[i]) >= threshold ? ceil(x[i]) : floor(x[i])
    end
end

# Hash function for cycle detection (only hashes integer variable values)
function hash_solution(
    x::Vector{Float64},
    integer_indices::Vector{Int}
)::Int

    hash(tuple((x[i] for i in integer_indices)...))
end

function perturb_solution!(
    x::Vector{Float64},
    x_round::Vector{Float64},
    binary::Vector{Int},
    integer::Vector{Int},
    lp_cols::Vector{Ptr{SCIP.SCIP_COL}},
)::Nothing

    for i in binary
        if rand() < DEF_PERTURB_FRACTION
            x_round[i] = 1.0 - x_round[i]
        end
    end

    for i in integer
        if rand() < DEF_PERTURB_FRACTION
            var = SCIP.SCIPcolGetVar(lp_cols[i])
            lb = SCIP.SCIPvarGetLbLocal(var)
            ub = SCIP.SCIPvarGetUbLocal(var)
            r = rand()

            newval = if (ub - lb) < DEF_BIGBIGM
                floor(lb + (1 + ub - lb) * r)
            elseif (x_round[i] - lb) < DEF_BIGM
                lb + (2 * DEF_BIGM - 1) * r
            elseif (ub - x_round[i]) < DEF_BIGM
                ub - (2 * DEF_BIGM - 1) * r
            else
                x[i] + (2 * DEF_BIGM - 1) * r - DEF_BIGM
            end

            x_round[i] = clamp(floor(newval), lb, ub)
        end
    end
end

# Helper function to check constraint feasibility
function check_feasibility(
    scip::Ptr{SCIP.SCIP_},
    solution::Vector{Float64},
    col_to_idx::Dict{Ptr{SCIP.SCIP_COL}, Int},
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)

    ptr_rows = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, ncols)

    # Check bounds
    for j in 1:ncols
        var = SCIP.SCIPcolGetVar(lp_cols[j])
        lb = SCIP.SCIPvarGetLbLocal(var)
        ub = SCIP.SCIPvarGetUbLocal(var)

        if solution[j] < lb - tolerance || solution[j] > ub + tolerance
            return false
        end
    end

    # Constraint check using rows
    for i in 1:nrows
        row = lp_rows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonz_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonz_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        activity = 0.0
        
        for k in 1:nnonz
            col = nonz_cols[k]
            idx = col_to_idx[col]
            activity += nonz_vals[k] * solution[idx]
        end
    
        constant = SCIP.SCIProwGetConstant(row)
        lhs = SCIP.SCIProwGetLhs(row) - constant
        rhs = SCIP.SCIProwGetRhs(row) - constant

        if lhs > -SCIP.SCIPinfinity(scip) && activity < lhs - tolerance
            return false
        end

        if rhs < SCIP.SCIPinfinity(scip) && activity > rhs + tolerance
            return false
        end
    end

    return true
end

# Helper function to check integrality
function check_integrality(
    solution::Vector{Float64}, 
    integer_indices::Vector{Int}, 
    tolerance::Float64=DEF_TOLERANCE
)::Bool

    for i in integer_indices
        if abs(solution[i] - round(solution[i])) > tolerance
            return false
        end
    end
    return true
end

function submit_solution_to_scip(
    scip::Ptr{SCIP.SCIP_},
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
    lp_cols::Vector{Ptr{SCIP.SCIP_COL}},
    solution::Vector{Float64},
    nvars::Int32
)::Bool

    sol_ptr = Ref{Ptr{SCIP.SCIP_SOL}}()
    SCIP.SCIPcreateSol(scip, sol_ptr, heur_ptr)
    sol = sol_ptr[]

    # Set solution values
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)
        SCIP.SCIPsetSolVal(scip, sol, var, solution[j])
    end

    # Try to add solution
    stored = Ref{SCIP.SCIP_Bool}()
    SCIP.SCIPtrySol(scip, sol, SCIP.TRUE, SCIP.FALSE, SCIP.TRUE, SCIP.TRUE, SCIP.TRUE, stored)

    if stored[] == SCIP.TRUE
        return true
    else
        SCIP.SCIPfreeSol(scip, sol_ptr)
        return false
    end
end