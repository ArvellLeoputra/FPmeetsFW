# Fix integer variables to xRound and solve the LP to adjust continuous variables
# Fallback when xRound alone is not accepted as a feasible MIP solution.
# `remainingTime` (seconds) bounds the dive LP: SCIPsolveDiveLP has no time argument, so we
# temporarily tighten `limits/time` (the diving solve derives its deadline from it) and restore
# it afterwards. A timed-out dive returns solstat != OPTIMAL, handled by the non-OPTIMAL branch.
function diveSolve(
    scip::Ptr{SCIP.SCIP_},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    intIdx::Vector{Int},
    xRound::Vector{Float64},
    ncols::Int32,
    remainingTime::Float64
)::Tuple{Bool, Vector{Float64}}
    # temporarily save the current SCIP time limit
    prevTimeLimit = Ref{SCIP.SCIP_Real}(0.0)
    SCIP.@SCIP_CALL SCIP.SCIPgetRealParam(scip, "limits/time", prevTimeLimit)

    SCIP.SCIPstartDive(scip)
    try
        SCIP.@SCIP_CALL SCIP.SCIPsetRealParam(scip, "limits/time", SCIP.SCIPgetSolvingTime(scip) + max(1e-3, remainingTime))

        # Fix integer variables to their rounded values
        for i in intIdx
            var = SCIP.SCIPcolGetVar(lpCols[i])
            SCIP.SCIPchgVarLbDive(scip, var, xRound[i])
            SCIP.SCIPchgVarUbDive(scip, var, xRound[i])
        end

        lperror = Ref{SCIP.SCIP_Bool}(SCIP.FALSE)
        cutoff  = Ref{SCIP.SCIP_Bool}(SCIP.FALSE)

        SCIP.SCIPsolveDiveLP(scip, -1, lperror, cutoff)

        if lperror[] == SCIP.TRUE || cutoff[] == SCIP.TRUE
            return false, Float64[]
        end

        if SCIP.SCIPgetLPSolstat(scip) == SCIP.SCIP_LPSOLSTAT_OPTIMAL
            sol = [SCIP.SCIPcolGetPrimsol(lpCols[j]) for j in 1:ncols]
            return true, sol
        else
            return false, Float64[]
        end
    finally
        SCIP.SCIPendDive(scip)
        # restore the SCIP time limit
        SCIP.@SCIP_CALL SCIP.SCIPsetRealParam(scip, "limits/time", prevTimeLimit[])
    end
end

# Submit a solution to SCIP, returning true if accepted, false otherwise
function submitSolution(
    scip::Ptr{SCIP.SCIP_},
    heurPtr::Ptr{SCIP.SCIP_HEUR},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    sol::Vector{Float64},
    ncols::Int32
)
    solPtr = Ref{Ptr{SCIP.SCIP_SOL}}()
    SCIP.SCIPcreateSol(scip, solPtr, heurPtr)
    scipSol = solPtr[]

    for j in 1:ncols
        col = lpCols[j]
        var = SCIP.SCIPcolGetVar(col)
        SCIP.SCIPsetSolVal(scip, scipSol, var, sol[j])
    end

    stored = Ref{SCIP.SCIP_Bool}()
    SCIP.SCIPtrySol(scip, scipSol, SCIP.FALSE, SCIP.FALSE, SCIP.TRUE, SCIP.TRUE, SCIP.TRUE, stored)

    if stored[] == SCIP.TRUE
        return true
    else
        SCIP.SCIPfreeSol(scip, solPtr)
        return false
    end
end