using SCIP
using JuMP
using FrankWolfe
using GLPK
import MathOptInterface


const MOI = MathOptInterface

mutable struct FPFWHeuristic <: SCIP.Heuristic
    called::Int64
    lmo::Union{Nothing, FrankWolfe.MathOptLMO}
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

    lb = zeros(SCIP.SCIP_Real, nvars)
    ub = zeros(SCIP.SCIP_Real, nvars)
    binary = []
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
    
    lmo = heur.lmo

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

SCIP.SCIPreadProb(scip, "filename.mps", C_NULL)  # Insert MPS file name here

# Build MOI model from SCIP model
nvars = SCIP.SCIPgetNVars(scip)
all_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetVars(scip), nvars)

moi_model = MOI.Utilities.Model{Float64}()
x = add_variables(moi_model, nvars)

# Add bounds
for j in 1:nvars
    var = all_vars[j]
    lb = SCIP.SCIPvarGetLbOriginal(var)
    ub = SCIP.SCIPvarGetUbOriginal(var)
    
    if lb > -SCIP.SCIPinfinity(scip)
        add_constraint(moi_model, x[j], GreaterThan(lb))
    end
    if ub < SCIP.SCIPinfinity(scip)
        add_constraint(moi_model, x[j], LessThan(ub))
    end
end

# Add original constraints
nconss = SCIP.SCIPgetNOrigConss(scip)
orig_conss = unsafe_wrap(Vector{Ptr{SCIP.SCIP_CONS}}, SCIP.SCIPgetOrigConss(scip), nconss)

for i in 1:nconss
    cons = orig_conss[i]
    conshdlr = SCIP.SCIPconsGetHdlr(cons)
    constype = unsafe_string(SCIP.SCIPconshdlrGetName(conshdlr))
    
    if constype == "linear"
        nvars_in_cons = SCIP.SCIPgetNVarsLinear(scip, cons)
        cons_vars = unsafe_wrap(Vector{Ptr{SCIP.SCIP_VAR}}, SCIP.SCIPgetVarsLinear(scip, cons), nvars_in_cons)
        cons_vals = unsafe_wrap(Vector{SCIP.SCIP_Real}, SCIP.SCIPgetValsLinear(scip, cons), nvars_in_cons)
        
        # Map SCIP variables to MOI variable indices
        var_to_idx = Dict(all_vars[k] => k for k in 1:nvars)
        
        terms = [ScalarAffineTerm(cons_vals[j], x[var_to_idx[cons_vars[j]]]) for j in 1:nvars_in_cons]
        aff = ScalarAffineFunction(terms, 0.0)
        
        lhs = SCIP.SCIPgetLhsLinear(scip, cons)
        rhs = SCIP.SCIPgetRhsLinear(scip, cons)
        
        if lhs > -SCIP.SCIPinfinity(scip)
            add_constraint(moi_model, aff, GreaterThan(lhs))
        end
        if rhs < SCIP.SCIPinfinity(scip)
            add_constraint(moi_model, aff, LessThan(rhs))
        end
    end
end

opt_model = GLPK.Optimizer()
copy_to(opt_model, moi_model)
lmo = FrankWolfe.MathOptLMO(opt_model)

heur = FPFWHeuristic(0, lmo)
SCIP.include_heuristic(
    backend, 
    heur,
    name="FPFWHeuristic", 
    priority=9999, 
    timing_mask=SCIP.SCIP_HEURTIMING_DURINGLPLOOP
)

SCIP.set_parameter(scip, "limits/nodes", 1)
SCIP.SCIPsolve(scip)