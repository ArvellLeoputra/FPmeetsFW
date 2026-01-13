include("dependencies.jl")
include("helper.jl")
include("lmo_builder.jl")

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

    println("\n" * "="^80)
    println("FPFW Heuristic Called (Attempt #$(heur.called))")
    println("="^80)

    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    println("There are $nvars variables and $nrows rows in the current LP.")

    # Build LMO from current LP
    if heur.lmo === nothing
        println("Building LMO...")
        println(typeof(scip))
        heur.lmo = build_lmo_from_scip_lp(scip, nvars, nrows)
        println("LMO built successfully")
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

    # println("LP solution: ", current_solution)
    
    for iter in 1:DEF_FP_MAX_ITER
        println("\n--- Iteration $iter ---")

        # Step 1: Rounding LP feasible solution
        x_round = copy(x)
        for i in binary
            x_round[i] = round(x[i])
        end
        
        # Step 2: Euclidean Projection using Frank-Wolfe
        x_new, _ = frank_wolfe(
            x -> f_proj(x, x_round),
            (g, x) -> g_proj(g, x, x_round),
            lmo,
            x,
            max_iteration = DEF_FW_MAX_ITER,
            verbose=false,
            line_search = FrankWolfe.Adaptive()
        )

        # Step 3: Check feasibility and integrality
        # Check integrality
        is_integral = check_integrality(x_new, binary)
        println("Integrality check: $is_integral")

        # Remove later -- just for debugging
        is_feasible = check_feasibility(scip, x_new)
        println("Feasibility check: $is_feasible")

        # Just for debugging
        if !is_feasible
            println("FrankWolfe return infeasible solution!")
        end

        if is_integral
            println("Found feasible integer solution after $iter iterations")
            
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
            distance = sqrt(sum((x_new[i] - x[i])^2 for i in binary))
            println("Distance moved: $(round(distance, digits=6))")

            if distance < DEF_TOLERANCE
                println("Converged to fixed point without finding integer solution")
                break
            end

            # If not feasible, continue with the non-rounded solution for next FW iteration
            println("Solution not feasible, continuing FW...")
            x .= x_new
        end
    end

    if result == SCIP.SCIP_DIDNOTFIND
        println("\nFailed to find feasible solution after $DEF_FP_MAX_ITER iterations")
    end

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