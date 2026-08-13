# Filters a directory of MIPLIB 2017 benchmark instances against the criteria used 
# to justify the benchmark selection in Section 5.3.
#
# For each instance this script solves ONLY the root node (node_limit=1) with
# all heuristics, cuts, and separators disabled (minimal_setup), with presolve
# enabled (matching the real FPFW run config in settings/fpfw.cfg), and flags
# instances that should be EXCLUDED from the FPFW benchmark if the root solve
# shows:
#   - no binary/general-integer variables (pure LP)
#   - root node already solved to optimality (nothing for a primal heuristic to do)
#   - infeasible or unbounded
#   - root LP/presolve did not finish within rootTimeLimit
#
# This intentionally reuses src/scip/setup.jl (minimal_setup) so the "no
# heuristics, no cuts, no separators" configuration matches the rest of the
# codebase exactly, rather than redefining it here.
#
# Usage (from the repo root):
#   julia --project scripts/filter_instances.jl <instance_dir> [output.csv] [rootTimeLimit] [mode]
#   mode: run (default, sequential + resumable) | task
#
# Example:
#   julia --project scripts/filter_instances.jl ./miplib_full filter_report.csv 300.0

using JuMP
using SCIP

include("../src/scip/setup.jl")

struct RootSolve
    readOk::Bool
    nVars::Int32
    nBin::Int32
    nInt::Int32
    nCont::Int32
    status::Symbol       # :optimal, :infeasible, :unbounded, :nodelimit, :timelimit, :unknown
    rootTime::Float64
end

struct InstanceReport
    name::String
    nVars::Int32
    nBin::Int32
    nInt::Int32
    nCont::Int32
    status::Symbol
    rootTime::Float64
    exclude::Bool
    reason::String
end

function classifyStatus(scip::SCIP.SCIPData)
    status = SCIP.SCIPgetStatus(scip)
    if status == SCIP.SCIP_STATUS_OPTIMAL
        return :optimal
    elseif status == SCIP.SCIP_STATUS_INFEASIBLE
        return :infeasible
    elseif status == SCIP.SCIP_STATUS_UNBOUNDED
        return :unbounded
    elseif status == SCIP.SCIP_STATUS_NODELIMIT
        return :nodelimit   # expected outcome: root node processed, branching would continue
    elseif status == SCIP.SCIP_STATUS_TIMELIMIT
        return :timelimit
    else
        return :unknown
    end
end

# Reads filePath fresh and solves ONLY the root node (node_limit=1) with all
# heuristics, cuts, and separators disabled (minimal_setup), either with or
# without presolve.
function solveRoot(filePath::String, rootTimeLimit::Float64, presolve::Bool)::RootSolve
    model = minimal_setup(time_limit=rootTimeLimit, node_limit=1, presolve=presolve)
    backend = JuMP.unsafe_backend(model)
    scip = backend.inner

    readOk = true
    try
        SCIP.SCIPreadProb(scip, filePath, C_NULL)
    catch
        readOk = false
    end

    if !readOk
        return RootSolve(false, Int32(0), Int32(0), Int32(0), Int32(0), :unknown, 0.0)
    end

    nVars = SCIP.SCIPgetNOrigVars(scip)
    nBin  = SCIP.SCIPgetNBinVars(scip)
    nInt  = SCIP.SCIPgetNIntVars(scip)
    nCont = SCIP.SCIPgetNContVars(scip)

    startT = time()
    SCIP.SCIPsolve(scip)
    rootTime = time() - startT

    status = classifyStatus(scip)

    return RootSolve(true, nVars, nBin, nInt, nCont, status, rootTime)
end

function checkInstance(filePath::String, rootTimeLimit::Float64)::InstanceReport
    name = basename(filePath)

    r = solveRoot(filePath, rootTimeLimit, true)

    if !r.readOk
        return InstanceReport(name, Int32(0), Int32(0), Int32(0), Int32(0),
            :unknown, 0.0, true, "failed to read")
    end

    exclude = false
    reason = ""

    if r.nBin + r.nInt == 0
        exclude = true
        reason = "no binary/general-integer variables"
    elseif r.status == :optimal
        exclude = true
        reason = "solved at root node without a primal heuristic"
    elseif r.status == :infeasible
        exclude = true
        reason = "infeasible"
    elseif r.status == :unbounded
        exclude = true
        reason = "unbounded"
    elseif r.status == :timelimit
        exclude = true
        reason = "root LP/presolve did not finish within rootTimeLimit=$(rootTimeLimit)s"
    end

    return InstanceReport(name, r.nVars, r.nBin, r.nInt, r.nCont, r.status, r.rootTime, exclude, reason)
end

function findInstanceFiles(instanceDir::String)
    files = filter(readdir(instanceDir)) do f
        endswith(f, ".mps") || endswith(f, ".lp") || endswith(f, ".mps.gz") || endswith(f, ".lp.gz")
    end
    sort!(files)
    return files
end

# Runs every instance in instanceDir sequentially (root node, presolve on) and writes the
# report incrementally (flushed after every instance, so progress is visible live and
# survives a crash/kill partway through).
#
# Resumable: if outCsv already has rows (from a previous killed/crashed run), already-processed
# instances are skipped and new results are appended.
function run(instanceDir::String, outCsv::String, rootTimeLimit::Float64)
    files = findInstanceFiles(instanceDir)

    done = alreadyProcessed(outCsv)
    remaining = filter(f -> !(f in done), files)
    fileMode = isempty(done) ? "w" : "a"

    if !isempty(done)
        println("Resuming: $(length(done))/$(length(files)) already done, $(length(remaining)) remaining.")
    end

    open(outCsv, fileMode) do io
        if fileMode == "w"
            println(io, "name,nVars,nBin,nInt,nCont,status,rootTime,exclude,reason")
            flush(io)
        end

        for (i, f) in enumerate(remaining)
            k = length(done) + i
            path = joinpath(instanceDir, f)
            print("[$k/$(length(files))] $f ... ")
            flush(stdout)

            r = checkInstance(path, rootTimeLimit)

            println(r.exclude ? "EXCLUDE ($(r.reason))" : "keep")
            flush(stdout)

            println(io, "$(r.name),$(r.nVars),$(r.nBin),$(r.nInt),$(r.nCont),$(r.status),$(round(r.rootTime,digits=3)),$(r.exclude),\"$(r.reason)\"")
            flush(io)
        end
    end

    dataLines = collect(eachline(outCsv))[2:end]
    totalExcluded = count(l -> !isempty(l) && split(l, ',')[8] == "true", dataLines)
    println("\nTotal instances checked: $(length(files))")
    println("Excluded: $totalExcluded")
    println("Kept: $(length(files) - totalExcluded)")
    println("Report written to $outCsv")
end

# Reads the name column of an already-written report CSV (if any), so a run
# can be resumed instead of restarted from scratch after a crash/kill.
function alreadyProcessed(outCsv::String)::Set{String}
    done = Set{String}()
    if isfile(outCsv)
        for (i, line) in enumerate(eachline(outCsv))
            i == 1 && continue   # header
            isempty(line) && continue
            name = split(line, ',')[1]
            push!(done, name)
        end
    end
    return done
end

# Runs a SINGLE instance (root node, presolve on, see checkInstance) selected by its
# 1-based position in the sorted instance directory listing, and writes exactly one
# CSV row (with header) to outDir/task_<taskIndex>.csv. Meant for SLURM array jobs:
# one array task = one instance, so all instances run in parallel instead of
# sequentially. Combine the per-task files afterward with scripts/merge_filter_tasks.sh.
function runTask(instanceDir::String, outDir::String, rootTimeLimit::Float64, taskIndex::Int)
    files = findInstanceFiles(instanceDir)
    if taskIndex < 1 || taskIndex > length(files)
        error("taskIndex=$taskIndex out of range 1:$(length(files))")
    end
    f = files[taskIndex]
    path = joinpath(instanceDir, f)

    println("[task $taskIndex/$(length(files))] $f ...")
    flush(stdout)

    r = checkInstance(path, rootTimeLimit)

    mkpath(outDir)
    outPath = joinpath(outDir, "task_$(taskIndex).csv")
    open(outPath, "w") do io
        println(io, "name,nVars,nBin,nInt,nCont,status,rootTime,exclude,reason")
        println(io, "$(r.name),$(r.nVars),$(r.nBin),$(r.nInt),$(r.nCont),$(r.status),$(round(r.rootTime,digits=3)),$(r.exclude),\"$(r.reason)\"")
    end

    println(r.exclude ? "EXCLUDE ($(r.reason))" : "keep")
    println("Task report written to $outPath")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        error("Usage: julia --project scripts/filter_instances.jl <instance_dir> [output.csv] [rootTimeLimit] [mode]\n" *
              "  mode: run (default, sequential + resumable) | task\n" *
              "  task mode: julia --project scripts/filter_instances.jl <instance_dir> <outDir> <rootTimeLimit> task <taskIndex>")
    end
    instanceDir = ARGS[1]
    outArg = length(ARGS) >= 2 ? ARGS[2] : "filter_report.csv"
    rootTimeLimit = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 300.0
    mode = length(ARGS) >= 4 ? ARGS[4] : "run"

    if mode == "run"
        run(instanceDir, outArg, rootTimeLimit)
    elseif mode == "task"
        length(ARGS) >= 5 || error("task mode requires a 5th argument: taskIndex")
        runTask(instanceDir, outArg, rootTimeLimit, parse(Int, ARGS[5]))
    else
        error("Unknown mode: $mode. Choose from run, task")
    end
end