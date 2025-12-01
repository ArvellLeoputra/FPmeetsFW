using SCIP
using JuMP
# using Printf
using FrankWolfe

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
end

struct MyLMO <: FrankWolfe.LinearMinimizationOracle
    model::Model
    x::Vector{VariableRef}
end

function MyLMO(A, b, l, u)
    m, n = size(A)

    model = Model(SCIP.Optimizer) # create inner model

    set_silent(model) # turn off output for the inner model

    @variable(model, l[j] <= x[j=1:n] <= u[j])

    # Linear equalities Ax >= b
    for i in 1:m
        @constraint(model, sum(A[i,j] * x[j] for j in 1:n) >= b[i])
    end

    @objective(model, Min, 0) # dummy objective, since JuMP requires one
    # The real objective will be set in compute_extreme_point

    return MyLMO(model, x)
end

function FrankWolfe.compute_extreme_point(lmo::MyLMO, direction; v=nothing, kwargs...)
    model = lmo.model
    x = lmo.x

    # Replace objective without rebuilding model
    @objective(model, Min, sum(direction[j] * x[j] for j in eachindex(x)))

    optimize!(model)
    println("LMO Called solution status: ", termination_status(model))
    return value.(x)
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

    # Get LP columns and rows
    ptr_cols = SCIP.SCIPgetLPCols(scip) # SCIP_COL**  // an array of pointers
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

    ptr_rows = SCIP.SCIPgetLPRows(scip) # SCIP_ROW**  // an array of pointers
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    cols_dict = Dict(lp_cols[k] => k for k in 1:nvars)  # Map each column to its index

    A = zeros(SCIP.SCIP_Real, nrows, nvars)
    b = zeros(SCIP.SCIP_Real, nrows)

    # Get the constraint matrix
    for i in 1:nrows
        row = lp_rows[i]

        constant = SCIP.SCIProwGetConstant(row)
        
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzero_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz)
        nonzero_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz)

        tmp_row = zeros(SCIP.SCIP_Real, nvars)
        for j in 1:nnonz
            k = cols_dict[nonzero_cols[j]]
            tmp_row[k] = nonzero_vals[j]
        end

        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

        if lhs != -SCIP.SCIPinfinity(scip)
            A[i, :] = tmp_row
            b[i] = lhs - constant

        elseif rhs != SCIP.SCIPinfinity(scip)
            A[i, :] = -tmp_row
            b[i] = constant - rhs

        else
            error("Row $i has both LHS and RHS infinite — cannot normalize to ≥ form")
        end
    end

    lb = zeros(SCIP.SCIP_Real, nvars)
    ub = zeros(SCIP.SCIP_Real, nvars)
    binary = []
    current_solution = zeros(SCIP.SCIP_Real, nvars)

    # Identify binary variables and get current LP solution
    for j in 1:nvars
        col = lp_cols[j]
        var = SCIP.SCIPcolGetVar(col)

        ub[j] = SCIP.SCIPcolGetUb(col)
        lb[j] = SCIP.SCIPcolGetLb(col)

        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binary, j)
        end

        current_solution[j] = SCIP.SCIPcolGetPrimsol(col) # Get current LP solution
    end

    # Penalty function with Manhattan norm
    function f(p)
        sum = 0.0
        for i in 1:length(p)
            if i in binary
                sum += min(p[i], 1.0 - p[i])
            end
        end
        return sum
    end

    function g(storage, p)
        storage .= 0.0
        for i in 1:length(p)
            if i in binary
                if p[i] < 0.5
                    storage[i] = 1.0
                else
                    storage[i] = -1.0
                end
            end
        end
        return storage
    end

    # Penalty function with Euclidean norm
    function f2(p)
        sum = 0.0
        for i in 1:length(p)
            if i in binary
                sum += min(p[i]^2, (1.0 - p[i])^2)
            end
        end
        return sum
    end

    function g2(storage, p)
        storage .= 0.0
        for i in 1:length(p)
            if i in binary
                if p[i] < 0.5
                    storage[i] = 2.0 * p[i]
                else
                    storage[i] = -2.0 * (1.0 - p[i])
                end
            end
        end
        return storage
    end

    lmo = MyLMO(A, b, lb, ub)
    
    p_opt, _ = frank_wolfe(
        f2,
        g2,
        lmo,
        current_solution,
        verbose=true,
        line_search = FrankWolfe.FixedStep(1.0)
    )
    return (SCIP.SCIP_OKAY, result)
end

model = direct_model(SCIP.Optimizer())
backend =  JuMP.unsafe_backend(model)
scip = backend.inner
heur = FPFWHeuristic(0)
SCIP.include_heuristic(
    backend, 
    heur,
    name="FPFWHeuristic", 
    priority=9999, 
    timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
)

SCIP.SCIPreadProb(scip, "filename.mps", C_NULL)  # Insert MPS file name here
SCIP.set_parameter(scip, "limits/nodes", 1)
SCIP.SCIPsolve(scip)
