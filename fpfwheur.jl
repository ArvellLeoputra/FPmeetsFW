using SCIP
using JuMP
using FrankWolfe

mutable struct FPFWHeuristic <: SCIP.Heuristic
    max_iter::Int
end

function SCIP.find_primal_solution(
    scip::Ptr{SCIP.SCIP_},
    heur::FPFWHeuristic,
    heurtiming::SCIP.SCIP_HEURTIMING,
    nodeinfeasible::SCIP.Bool,
    heur_ptr::Ptr{SCIP.SCIP_HEUR},
)::Tuple{SCIP.SCIP_RETCODE, SCIP.SCIP_RESULT}

    @info("FWFP heuristic called")
    @assert SCIP.SCIPhasCurrentNodeLP(scip) == SCIP.TRUE  # always true, since we set timing to DURINGLPLOOP
    result = SCIP.SCIP_DIDNOTRUN

    n_vars = SCIP.SCIPgetNLPCols(scip)
    n_rows = SCIP.SCIPgetNLPRows(scip)
    @info "There are $n_vars variables and $n_rows rows in the current LP."

    # For each column, there is a single SCIP_COL object, while for each row there is a single SCIP_ROW object.
    # First, we get all the SCIP_ROW objects
    # We have an object called SCIP_ROW, this lies in the heap
    # So to keep track of the object we have a pointer to the object SCIP_ROW*
    # Now we an array of such pointers, one for each row so we have SCIP_ROW** (pointer to pointer)
    ptr_rows = SCIP.SCIPgetLPRows(scip) # SCIP_ROW** // an array of pointers
    row_pointers = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, n_rows)
    # row_pointers is an array of pointers to SCIP_ROW objects

    # gets all constraints printed out
    for (i, row_ptr) in enumerate(row_pointers)
        # gets the LHS and RHS of the row, i.e., the constraint bounds
        rhs = SCIP.SCIProwGetRhs(row_ptr)
        lhs = SCIP.SCIProwGetLhs(row_ptr)

        shift = SCIP.SCIProwGetConstant(row_ptr)
        nnonz = SCIP.SCIProwGetNNonz(row_ptr)

        non_zero_cols = SCIP.SCIProwGetCols(row_ptr) # SCIP_COL**
        non_zero_cols_ptrs = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, non_zero_cols, nnonz)
        non_zero_vals = SCIP.SCIProwGetVals(row_ptr) # double*
        non_zero_vec = unsafe_wrap(Vector{SCIP.SCIP_Real}, non_zero_vals, nnonz)

        print("$lhs <=")
        for j in 1:nnonz
            var = SCIP.SCIPcolGetVar(non_zero_cols_ptrs[j])
            name = unsafe_string(SCIP.SCIPvarGetName(var))
            coeff = non_zero_vec[j]
            print(" $coeff*$name ")
        end
        print(" + $shift <=  $rhs\n")
    end

    return (SCIP.SCIP_OKAY, SCIP.SCIP_DIDNOTFIND)
end

optimizer = SCIP.Optimizer()

model = direct_model(optimizer)
set_attribute(model, "presolving/maxrounds", 0)

@variable(model, x[1:3] >= 0)
@variable(model, 3 >= y >= 2, Int)
@objective(model, Max, x[1] + 2 * x[2] + 3 * x[3] + y)
@constraint(model, -x[1] + x[2] + x[3] + 10 * y <= 20)
@constraint(model, x[1] - 3 * x[2] + x[3] <= 30)
@constraint(model, x[2] - 3.5 * y == 0)
@constraint(model, x[1] <= 40)
@constraint(model, 2 <= y <= 3)

scip = JuMP.unsafe_backend(model).inner
SCIP.include_heuristic(
    optimizer, 
    FPFWHeuristic(100);
    name="FPFWHeuristic", 
    description="Frank–Wolfe Feasibility Pump heuristic", 
    priority=100000, 
    frequency=0, 
    timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
)

optimize!(model)
assert_is_solved_and_feasible(model)
solution_summary(model)