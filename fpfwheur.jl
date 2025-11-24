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

function MyLMO(A, l, u)
    m, n = size(A)

    model = Model(SCIP.Optimizer)

    set_silent(model)

    # Variables with bounds
    @variable(model, l[j] <= x[j=1:n] <= u[j])

    # Linear equalities A*x = 0
    for i in 1:m
        @constraint(model, sum(A[i,j] * x[j] for j in 1:n) == 0)
    end

    # Dummy objective (will overwrite later)
    @objective(model, Min, 0)

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

    @info("FPFW heuristic called")
    nvars = SCIP.SCIPgetNLPCols(scip)
    nrows = SCIP.SCIPgetNLPRows(scip)
    @info "There are $nvars variables and $nrows rows in the current LP."
    constraint_matrix = zeros(SCIP.SCIP_Real, nrows, nvars + nrows)

    ptr_cols = SCIP.SCIPgetLPCols(scip) # SCIP_COL**  // an array of pointers
    lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, ptr_cols, nvars)

    ptr_rows = SCIP.SCIPgetLPRows(scip) # SCIP_ROW**  // an array of pointers
    lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, ptr_rows, nrows)

    # Get the constraint matrix
    for i in 1:nrows
        row = lp_rows[i]
        nnonz = SCIP.SCIProwGetNNonz(row)
        nonzero_cols = unsafe_wrap(
            Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz
            )
        nonzero_vals = unsafe_wrap(
            Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz
        )

        for j in 1:nnonz
            for k in 1:nvars
                if lp_cols[k] == nonzero_cols[j]
                    constraint_matrix[i, k] = nonzero_vals[j]
                end
            end
        end

        constraint_matrix[i, nvars + i] = 1  # Slack Variable
    end

    binary = []
    solution = zeros(SCIP.SCIP_Real, nvars + nrows)

    for i in 1:nrows
        row = lp_rows[i]
        if SCIP.SCIProwGetConstant(row) != 0
            @error "Row with nonzero constant"
        end
    end

    for i in 1:nvars
        col = lp_cols[i]
        var = SCIP.SCIPcolGetVar(col)
        if SCIP.SCIPvarIsBinary(var) == SCIP.TRUE
            push!(binary, i)
        end
        solution[i] = SCIP.SCIPcolGetPrimsol(col)
    end

    upper_bounds = zeros(SCIP.SCIP_Real, nvars + nrows)
    lower_bounds = zeros(SCIP.SCIP_Real, nvars + nrows)

    for i in 1:nvars
        col = lp_cols[i]
        upper_bounds[i] = SCIP.SCIPcolGetUb(col)
        lower_bounds[i] = SCIP.SCIPcolGetLb(col)
    end

    for i in 1:nrows
        upper_bounds[nvars + i] = -SCIP.SCIProwGetLhs(lp_rows[i])
        lower_bounds[nvars + i] = -SCIP.SCIProwGetRhs(lp_rows[i])
        # SCIP.SCIPprintRow(scip, lp_rows[i], C_NULL)
    end

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

    gradient_storage = zeros(SCIP.SCIP_Real, nvars + nrows)
    lmo = MyLMO(constraint_matrix, lower_bounds, upper_bounds)
    p_opt, _ = frank_wolfe(
        f,
        g,
        lmo,
        solution,
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

SCIP.SCIPreadProb(scip, "MPSFILE", C_NULL)  # Insert MPS file name here
SCIP.set_parameter(scip, "limits/nodes", 1)
SCIP.SCIPsolve(scip)
