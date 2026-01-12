using SCIP
using JuMP
using FrankWolfe
using GLPK
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

# Helper function to check constraint feasibility
function check_feasibility(
    scip::Ptr{SCIP.SCIP_}, 
    solution::Vector{Float64},
    tolerance=1e-6
)::Bool

    ncols = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    all_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetVars(scip), ncols)

    ptr_rows = SCIP.SCIPgetLPRows(scip)
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    ptr_cols = SCIP.SCIPgetLPCols(scip)
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, ncols)

    # Check bounds
    for j in 1:ncols
        vars = all_vars[j]
        lb = SCIP.SCIPvarGetLbLocal(vars)
        ub = SCIP.SCIPvarGetUbLocal(vars)

        if solution[j] < lb - tolerance || solution[j] > ub + tolerance
            # println("Variable $j violates bounds: $lb <= $(solution[j]) <= $ub")
            return false
        end
    end
    
    # Constraint check using rows
    for i in 1:nrows
        row = lp_rows[i]

        nnonz = SCIP.SCIProwGetNNonz(row)
        nonz_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonz_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        col_to_idx = Dict(lp_cols[k] => k for k in 1:ncols)
        
        activity = 0.0
        
        for k in 1:nnonz
            col = nonz_cols[k]
            idx = col_to_idx[col]
            activity += nonz_vals[k] * solution[idx]
        end
    
        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

        if lhs > -SCIP.SCIPinfinity(scip) && activity < lhs - tolerance
            # println("Constraint $i violates LHS: $activity >= $lhs")
            return false
        end

        if rhs < SCIP.SCIPinfinity(scip) && activity > rhs + tolerance
            #println("Constraint $i violates RHS: $activity <= $rhs")
            return false
        end 
    end

    return true
end

# Helper function to check integrality
function check_integrality(
    solution::Vector{Float64}, 
    binary_indices::Vector{Int}, 
    tolerance=1e-6
)::Bool

    for i in binary_indices
        if abs(solution[i] - round(solution[i])) > tolerance
            return false
        end
    end
    return true
end

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

    println("FPFW heuristic called")

    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    println("There are $nvars variables and $nrows rows in the current LP.")

    # Build LMO from current LP
    if heur.lmo == nothing
        println("Building LMO from current LP with $nvars variables...")

        ptr_cols = SCIP.SCIPgetLPCols(scip)
        lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

        ptr_rows = SCIP.SCIPgetLPRows(scip)
        lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

        # Build MOI model from current LP
        moi_model = MOI.Utilities.Model{Float64}()
        x = add_variables(moi_model, nvars)

        # Add variable domains from LP columns (lb <= x <= ub)
        for j in 1:nvars
            col = lp_cols[j]
            var = SCIP.SCIPcolGetVar(col)

            lb = SCIP.SCIPvarGetLbLocal(var)
            ub = SCIP.SCIPvarGetUbLocal(var)

            if lb > -SCIP.SCIPinfinity(scip)
                add_constraint(moi_model, x[j], GreaterThan(lb))
            end
            if ub < SCIP.SCIPinfinity(scip)
                add_constraint(moi_model, x[j], LessThan(ub))
            end
        end

        # Add row constraints (lhs <= Ax <= rhs)
        col_to_idx = Dict(lp_cols[k] => k for k in 1:nvars)
        
        for i in 1:nrows
            row = lp_rows[i]
            nnonz = SCIP.SCIProwGetNNonz(row)
            nonzero_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
            nonzero_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)
            
            terms = [ScalarAffineTerm(nonzero_vals[k], x[col_to_idx[nonzero_cols[k]]]) for k in 1:nnonz]
            aff = ScalarAffineFunction(terms, 0.0)
            
            lhs = SCIP.SCIProwGetLhs(row)
            rhs = SCIP.SCIProwGetRhs(row)
            
            if lhs > -SCIP.SCIPinfinity(scip)
                add_constraint(moi_model, aff, GreaterThan(lhs))
            end
            if rhs < SCIP.SCIPinfinity(scip)
                add_constraint(moi_model, aff, LessThan(rhs))
            end
        end

        opt_model = GLPK.Optimizer()
        copy_to(opt_model, moi_model)
        heur.lmo = FrankWolfe.MathOptLMO(opt_model)
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

    println("Number of variables: $nvars")
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

    # Penalty function with Euclidean norm
    function f2(p)
        sum = 0.0
        for i in binary
            sum += min(p[i]^2, (1.0 - p[i])^2)
        end
        return sum
    end

    function g2(storage, p)
        storage .= 0.0
        for i in binary
            if p[i] < 0.5
                storage[i] = 2.0 * p[i]
            else
                storage[i] = -2.0 * (1.0 - p[i])
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

    max_iterations = 100
    x = copy(current_solution)
    println("LP solution: ", current_solution)
    
    for iter in 1:max_iterations
        println("\n--- Iteration $iter ---")
        
        x_round = copy(x)
        for i in binary
            if x_round[i] < 0.5
                x_round[i] = 0.0
            else
                x_round[i] = 1.0
            end
        end
        
        x_fw, _ = frank_wolfe(
            x -> f_proj(x, x_round),
            (g, x) -> g_proj(g, x, x_round),
            lmo,
            x,
            max_iteration = 10,
            verbose=false,
            line_search = FrankWolfe.Adaptive()
        )

        println("Current FW solution: ", x_fw)

        # Check integrality
        is_integral = check_integrality(x_fw, binary)
        println("Integrality check: $is_integral")

        # Remove later -- just for debugging
        is_feasible = check_feasibility(scip, x_fw)
        println("Feasibility check: $is_feasible")

        if !is_integral
            distance = sqrt(sum((x_fw[i] - x[i])^2 for i in binary))
            println("Distance moved: $distance")
            if distance < 1e-6
                println("Converged to fixed point without finding integer solution")
                break
            end
        end

        if is_integral
            println("Found feasible integer solution after $iter iterations")
            
            # Add solution to SCIP
            sol_ptr = Ref{Ptr{SCIP.SCIP_SOL}}()
            SCIP.SCIPcreateSol(scip, sol_ptr, heur_ptr)
            sol = sol_ptr[]
            
            # Get all variables
            all_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetVars(scip), nvars)
            
            # Map LP columns back to variables
            for j in 1:nvars
                col = lp_cols[j]
                var = SCIP.SCIPcolGetVar(col)
                SCIP.SCIPsetSolVal(scip, sol, var, x_fw[j])
            end
            
            stored = Ref{SCIP.SCIP_Bool}()
            ret = SCIP.SCIPtrySol(scip, sol, SCIP.FALSE, SCIP.FALSE, SCIP.TRUE, SCIP.TRUE, SCIP.TRUE, stored)
            
            if stored[] == SCIP.TRUE
                obj = SCIP.SCIPgetSolOrigObj(scip, sol)
                println("Objective value of accepted solution: $obj")
                println("Solution accepted by SCIP!")
                result = SCIP.SCIP_FOUNDSOL
            else
                println("Solution rejected by SCIP")
                SCIP.SCIPfreeSol(scip, sol_ptr)
            end
            
            break
        else
            # If not feasible, continue with the non-rounded solution for next FW iteration
            println("Solution not feasible, continuing FW...")
            x .= x_fw
        end
    end

    if result == SCIP.SCIP_DIDNOTFIND
        println("\nFailed to find feasible solution after $max_iterations iterations")
    end

    return (SCIP.SCIP_OKAY, result)
end

model = direct_model(SCIP.Optimizer())
backend =  JuMP.unsafe_backend(model)
scip = backend.inner

SCIP.SCIPreadProb(scip, "filename.mps", C_NULL)  # Insert MPS file name here

heur = FPFWHeuristic(0, nothing)
SCIP.include_heuristic(
    backend, 
    heur,
    name="FPFWHeuristic", 
    priority=9999, 
    timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
)

SCIP.set_parameter(scip, "limits/nodes", 1)
SCIP.SCIPsolve(scip)