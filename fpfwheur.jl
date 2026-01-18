# Main FPFW Heuristic Implementation
function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}

    @assert SCIP.SCIPhasCurrentNodeLP(scip) == SCIP.TRUE  # always true, since we set timing to DURINGLPLOOP
    result = SCIP.SCIP_DIDNOTFIND

    heur.called += 1
    if heur.called > 1
        return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTRUN)
    end

    stats = FPFWStats()
    total_start_time = time()

    println("\n" * "="^80)
    println("FPFW Heuristic Called (Attempt #$(heur.called))")
    println("="^80)

    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    println("There are $nvars variables and $nrows rows in the current LP.")

    # Build LMO from current LP
    if heur.lmo === nothing
        lmo_start = time()
        println("Building LMO...")
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
        lmo_time = time() - lmo_start
        println("LMO built successfully in $(round(lmo_time, digits=3)) seconds.")
    end

    lmo = heur.lmo

    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

    binary = Int[]
    current_solution = zeros(SCIP.SCIP_Real, nvars)

    # Identify binary variables and get current LP solution
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)

        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binary, j)
        end

        current_solution[j] = SCIP.SCIPcolGetPrimsol(col) # Get current LP solution
    end

    println("Binary variables: $(length(binary)) out of $nvars")

    lp_obj = 0.0
    for j in 1:nvars
        lp_obj += current_solution[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j]))
    end

    println("Initial LP objective: $(round(lp_obj, digits=4))")

    # Penalty function with Manhattan norm
    function f(p)
        sum = 0.0
        for i in binary
            sum += min(p[i], 1.0 - p[i])
        end
        return sum
    end

    function g(storage, p)
        storage .= 0.0
        for i in binary
            if p[i] < 0.5
                storage[i] = 1.0
            else
                storage[i] = -1.0
            end
        end
        return storage
    end

    function f_proj(x, x_round)
        sum = 0.0
        for i in binary
            d = x[i] - x_round[i]
            sum += d * d
        end
        return 0.5 * sum
    end
    
    function g_proj(storage, x, x_round)
        storage .= 0.0
        for i in binary
            storage[i] = x[i] - x_round[i]
        end
        return storage
    end
    
    # Main Feasibility Pump with Frank-Wolfe Loop
    x = copy(current_solution)  # Initial LP feasible solution
    prev_x = copy(x)
    cycle = false

    # println("LP solution: ", current_solution)
    
    # Track best solution found
    best_solution = nothing
    best_objective = Inf

    for iter in 1:DEF_FP_MAX_ITER
        stats.fp_iterations = iter

        # Check time limit (5 minutes)
        elapsed = time() - total_start_time
        if elapsed > DEF_TIME_LIMIT
            println("Time limit reached ($(round(elapsed, digits=2))s > $(DEF_TIME_LIMIT)s)")
            break
        end

        # Step 1: Rounding LP feasible solution
        x_round = copy(x)
        for i in binary
            x_round[i] = round(x[i])
        end

        fw_traj = Vector{NTuple{5,Float64}}()
        fw_callback = FrankWolfe.make_print_callback(
            FrankWolfe.make_trajectory_callback(nothing, fw_traj),
            10,
            ("iter", "primal", "dual", "gap", "time"),
            "%6i %12e %12e %12e %12e\n",
            FrankWolfe.callback_state
        )

        stats.fw_calls += 1
        
        # Step 2: Euclidean Projection using Frank-Wolfe
        @time x_new, _ = frank_wolfe(
            y -> f_proj(y, x_round),
            (g, y) -> g_proj(g, y, x_round),
            lmo,
            x,
            max_iteration = DEF_FW_MAX_ITER,
            verbose = false,
            line_search = FrankWolfe.Adaptive(),
            epsilon = 1e-7,
            callback = fw_callback
        )

        fw_iters = length(fw_traj)
        stats.fw_iterations += fw_iters

        if fw_iters > 0
            stats.fw_time += fw_traj[end][5] - fw_traj[1][5]
        end

        # Step 3: Check feasibility and integrality and calculate current objective value
        is_integral = check_integrality(x_new, binary)
        is_feasible = check_feasibility(scip, x_new)
        distance = sqrt(sum((x_new[i] - prev_x[i])^2 for i in binary))

        obj_val = 0.0
        for j in 1:nvars
            obj_val += x_new[j] * SCIP.SCIPvarGetObj(SCIP.SCIPcolGetVar(lp_cols[j]))
        end

        # Just for debugging
        if !is_feasible
            println("FrankWolfe return infeasible solution!")
            break
        end

        if is_integral
            println("Found feasible integer solution at iteration $iter with objective $obj_val")
            stats.solution_found = true
            stats.final_objective = obj_val

            # Track best solution (for minimization problems)
            if obj_val < best_objective
                best_objective = obj_val
                best_solution = copy(x_new)
            end

            # Submit solution to SCIP
            if submit_solution_to_scip(scip, heur_ptr, lp_cols, x_new, nvars)
                result = SCIP.SCIP_FOUNDSOL
            else
                println("Solution rejected by SCIP")
                break
            end
        else
            # Short cycle detection
            # TODO: Implement cycle detection, maybe with a cache of previous rounded solutions
            if distance < DEF_TOLERANCE
                println("Converged to fixed point without finding integer solution")
                cycle = true
                break
            end

            # If not feasible, continue with the non-rounded solution for next FW iteration
            prev_x .= copy(x)
            x .= x_new
        end
    end

    # Print final summary
    total_time = time() - total_start_time
    println("\n" * "="^80)
    println("FPFW HEURISTIC SUMMARY")
    println("="^80)
    println("Total time:        $(round(total_time, digits=2)) seconds")
    println("FP iterations:     $(stats.fp_iterations)")
    println("FW calls:          $(stats.fw_calls)")
    println("FW iterations:     $(stats.fw_iterations)")
    println("Solution found:    $(stats.solution_found)")

    if best_solution !== nothing
        println("\nBEST SOLUTION FOUND:")
        println("  Objective value: $best_objective")
        println("  Solution vector (first 20 vars): $(best_solution[1:min(20, length(best_solution))])")
    else
        println("\nNo feasible integer solution found.")
    end
    println("="^80 * "\n")

    return (SCIP.SCIP_OKAY, result)
end

function submit_solution_to_scip(
    scip::Ptr{SCIP.SCIP_},
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
    lp_cols::Vector{Ptr{SCIP.SCIP_COL}},
    solution::Vector{Float64},
    nvars::Int
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
        println("Solution accepted!")
        return true
    else
        SCIP.SCIPfreeSol(scip, sol_ptr)
        return false
    end
end