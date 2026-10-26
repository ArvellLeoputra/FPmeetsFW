# Track per-variable round stability and lock/unlock the LMO's bounds as variables settle or change.
function updateLockState!(
    scip::Ptr{SCIP.SCIP_},
    st::StageState,
    config::FPFWConfig,
    lmo::FrankWolfe.LinearMinimizationOracle,
    xRound::Vector{Float64},
    prevRound::Vector{Float64},
    stats::FPFWStats,
    ncols::Int32
)
    if config.fwLockThreshold == 0
        return
    end

    # Skip the stage's first iteration: prevRound isn't a real prior observation yet, to avoid false positive
    if st.stageIter <= 1
        return
    end

    toUnlock = Int[]
    lockQueue = Int[]

    for i in st.activeIntIdx
        if SCIP.SCIPisEQ(scip, xRound[i], prevRound[i]) == SCIP.TRUE
            st.roundStableCount[i] += 1
            if st.roundStableCount[i] >= config.fwLockThreshold && i ∉ st.lockedIdx
                push!(lockQueue, i)
            end
        else
            st.roundStableCount[i] = 0
            if i in st.lockedIdx
                push!(toUnlock, i)
                delete!(st.lockedIdx, i)
            end
        end
    end

    if !isempty(toUnlock)
        LPIunlockVars!(lmo, toUnlock)
    end

    # Lock candidates one at a time, most-stable first, each feasibility-checked before being kept.
    # Capped per iteration so locking ramps up gradually.
    if !isempty(lockQueue)
        sort!(lockQueue, by = i -> st.roundStableCount[i], rev = true)
        cap = min(max(1, round(Int, config.fwLockMaxFrac * ncols)), config.fwLockCap)
        nNewLocks = 0
        for i in lockQueue
            if nNewLocks >= cap
                break
            end
            simplexItersBefore = lmo.simplexIters[]
            try
                LPIlockVars!(lmo, [i], [xRound[i]])
                if LPIcheckFeasible!(lmo)
                    push!(st.lockedIdx, i)
                    nNewLocks += 1
                else
                    LPIunlockVars!(lmo, [i])
                    stats.lockRejected += 1
                end
            catch e
                (e isa LMODeadlineExceeded || e isa LMOSolverError) || rethrow()
                # Out of time or a solver glitch mid-check: i is pinned but was never added to
                # lockedIdx, so unlockAll! would never find it later. Unlock it directly before giving up.
                LPIunlockVars!(lmo, [i])
                stats.lockCheckSimplexIters += lmo.simplexIters[] - simplexItersBefore
                break
            end
            stats.lockCheckSimplexIters += lmo.simplexIters[] - simplexItersBefore
        end
    end

    stats.lockVarsPeak = max(stats.lockVarsPeak, length(st.lockedIdx))
end

# Safety net for LMOPrimalInfeasible: undo all locks and reset the stability counters.
function unlockAll!(st::StageState, lmo::FrankWolfe.LinearMinimizationOracle)
    if isempty(st.lockedIdx)
        return
    end

    LPIunlockVars!(lmo, collect(st.lockedIdx))
    empty!(st.lockedIdx)
    st.roundStableCount .= 0
end
