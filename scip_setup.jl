# SCIP configuration and setup utilities for isolated testing of FPFW heuristic
function set_param(model::JuMP.Model, param::String, value)
    JuMP.set_attribute(model, param, value)
end

function disable_all_heuristics!(model::JuMP.Model)
    set_param(model, "heuristics/padm/freq", -1)
    set_param(model, "heuristics/ofins/freq", -1)
    set_param(model, "heuristics/trivialnegation/freq", -1)
    set_param(model, "heuristics/reoptsols/freq", -1)
    set_param(model, "heuristics/trivial/freq", -1)
    set_param(model, "heuristics/clique/freq", -1)
    set_param(model, "heuristics/locks/freq", -1)
    set_param(model, "heuristics/vbounds/freq", -1)
    set_param(model, "heuristics/shiftandpropagate/freq", -1)
    set_param(model, "heuristics/completesol/freq", -1)
    set_param(model, "heuristics/simplerounding/freq", -1)
    set_param(model, "heuristics/randrounding/freq", -1)
    set_param(model, "heuristics/zirounding/freq", -1)
    set_param(model, "heuristics/rounding/freq", -1)
    set_param(model, "heuristics/shifting/freq", -1)
    set_param(model, "heuristics/intshifting/freq", -1)
    set_param(model, "heuristics/oneopt/freq", -1)
    set_param(model, "heuristics/indicator/freq", -1)
    set_param(model, "heuristics/adaptivediving/freq", -1)
    set_param(model, "heuristics/farkasdiving/freq", -1)
    set_param(model, "heuristics/feaspump/freq", -1)
    set_param(model, "heuristics/conflictdiving/freq", -1)
    set_param(model, "heuristics/pscostdiving/freq", -1)
    set_param(model, "heuristics/fracdiving/freq", -1)
    set_param(model, "heuristics/nlpdiving/freq", -1)
    set_param(model, "heuristics/veclendiving/freq", -1)
    set_param(model, "heuristics/distributiondiving/freq", -1)
    set_param(model, "heuristics/objpscostdiving/freq", -1)
    set_param(model, "heuristics/rootsoldiving/freq", -1)
    set_param(model, "heuristics/linesearchdiving/freq", -1)
    set_param(model, "heuristics/guideddiving/freq", -1)
    set_param(model, "heuristics/rens/freq", -1)
    set_param(model, "heuristics/alns/freq", -1)
    set_param(model, "heuristics/rins/freq", -1)
    set_param(model, "heuristics/gins/freq", -1)
    set_param(model, "heuristics/lpface/freq", -1)
    set_param(model, "heuristics/crossover/freq", -1)
    set_param(model, "heuristics/undercover/freq", -1)
    set_param(model, "heuristics/subnlp/freq", -1)
    set_param(model, "heuristics/mpec/freq", -1)
    set_param(model, "heuristics/multistart/freq", -1)
    set_param(model, "heuristics/trysol/freq", -1)
end

function disable_separators!(model::JuMP.Model)
    set_param(model, "separating/disjunctive/freq", -1)
    set_param(model, "separating/impliedbounds/freq", -1)
    set_param(model, "separating/gomory/freq", -1)
    set_param(model, "separating/strongcg/freq", -1)
    set_param(model, "separating/aggregation/freq", -1)
    set_param(model, "separating/clique/freq", -1)
    set_param(model, "separating/zerohalf/freq", -1)
    set_param(model, "separating/mcf/freq", -1)
    set_param(model, "separating/flowcover/freq", -1)
    set_param(model, "separating/cmir/freq", -1)
    set_param(model, "separating/rapidlearning/freq", -1)
    set_param(model, "constraints/cardinality/sepafreq", -1)
    set_param(model, "constraints/SOS1/sepafreq", -1)
    set_param(model, "constraints/SOS2/sepafreq", -1)
    set_param(model, "constraints/varbound/sepafreq", -1)
    set_param(model, "constraints/knapsack/sepafreq", -1)
    set_param(model, "constraints/setppc/sepafreq", -1)
    set_param(model, "constraints/linking/sepafreq", -1)
    set_param(model, "constraints/or/sepafreq", -1)
    set_param(model, "constraints/and/sepafreq", -1)
    set_param(model, "constraints/xor/sepafreq", -1)
    set_param(model, "constraints/linear/sepafreq", -1)
    set_param(model, "constraints/orbisack/sepafreq", -1)
    set_param(model, "constraints/symresack/sepafreq", -1)
    set_param(model, "constraints/logicor/sepafreq", -1)
    set_param(model, "constraints/cumulative/sepafreq", -1)
    set_param(model, "constraints/nonlinear/sepafreq", -1)
    set_param(model, "separating/mixing/freq", -1)
    set_param(model, "separating/rlt/freq", -1)
    set_param(model, "constraints/indicator/sepafreq", -1)
end

function disable_cuts!(model::JuMP.Model)
    set_param(model, "separating/maxrounds", 0)
    set_param(model, "separating/maxroundsroot", 0)
end

function disable_presolving!(model::JuMP.Model)
    set_param(model, "presolving/maxrounds", 0)
end

function set_limits!(model::JuMP.Model; time_limit, node_limit)
    set_param(model, "limits/time", time_limit)
    set_param(model, "limits/nodes", node_limit)
end

function disable_root_node_propagation!(model::JuMP.Model)
    set_param(model, "propagating/maxroundsroot", 0)
end

function disable_strong_branching_lookahead!(model::JuMP.Model)
    set_param(model, "branching/relpscost/initcand", 0)
    set_param(model, "branching/relpscost/sbiterofs", 0)
    set_param(model, "branching/relpscost/sbiterquot", 0)
    set_param(model, "branching/relpscost/priority", -99999)
    set_param(model, "branching/random/priority", 100000)
end

function set_verbosity!(model::JuMP.Model, level::Int)
    # 0=none, 1=errors, 2=warnings, 3=normal, 4=high, 5=full
    set_param(model, "display/verblevel", level)
end

function minimal_setup(;
    time_limit=DEF_SCIP_TIME_LIMIT,
    node_limit=1,
    verbosity=0,
    presolve=false
)
    model = Model(SCIP.Optimizer)

    disable_all_heuristics!(model)
    disable_cuts!(model)
    disable_root_node_propagation!(model)

    if !presolve                                                                                                                                                                                                                  
        disable_presolving!(model)                                                                                                                                                                                              
    end

    set_verbosity!(model, verbosity)
    set_limits!(model, time_limit=time_limit, node_limit=node_limit)

    return model
end