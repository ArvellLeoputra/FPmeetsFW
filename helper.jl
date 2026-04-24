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
        "integer feasible solution accepted by SCIP at iteration $(stats.iter_found_solution)"
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
    println("FP iterations:     $(stats.fp_iterations)")
    println("FW iterations:     $(stats.fw_iterations)")
    println("FW time:           $(round(stats.fw_time, digits=2))s")
    println("Restarts:          $(stats.restarts)")
    println("Solution found:    $(stats.solution_found)")
    println("Exit reason:       $exit_msg")
    println("="^80 * "\n")
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
    
        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

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