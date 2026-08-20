# Combines N repeated filter_report_run<k>.csv files (same instances, run independently
# on the cluster N times) into one comparison table per instance, to see how much
# rootTime varies run-to-run before deciding on a final keep/exclude threshold.
#
# This does NOT make the final keep/exclude call itself -- it just lays out the raw
# numbers (all N rootTimes, mean, min, max, how many hit the 300s timelimit) so that
# decision can be made deliberately afterward.
#
# Usage:
#   julia --project scripts/aggregate_filter_runs.jl <outCsv> <run1.csv> <run2.csv> ... <runN.csv>
#
# Example:
#   julia --project scripts/aggregate_filter_runs.jl misc/filter_report_aggregated.csv \
#       misc/filter_report_run1.csv misc/filter_report_run2.csv misc/filter_report_run3.csv \
#       misc/filter_report_run4.csv misc/filter_report_run5.csv

struct Row
    nVars::String
    nBin::String
    nInt::String
    nCont::String
    status::String
    rootTime::Float64
end

function readReport(path::String)::Dict{String, Row}
    rows = Dict{String, Row}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue   # header
        isempty(line) && continue
        parts = split(line, ',')
        name = parts[1]
        rows[name] = Row(parts[2], parts[3], parts[4], parts[5], parts[6], parse(Float64, parts[7]))
    end
    return rows
end

function main()
    length(ARGS) >= 3 || error("Usage: julia --project scripts/aggregate_filter_runs.jl <outCsv> <run1.csv> <run2.csv> ...")

    outCsv = ARGS[1]
    runPaths = ARGS[2:end]
    n = length(runPaths)

    reports = [readReport(p) for p in runPaths]

    allNames = sort(collect(keys(reports[1])))
    for r in reports[2:end]
        allNames == sort(collect(keys(r))) || @warn "Instance sets differ between run reports -- some names may be missing in the output"
    end

    open(outCsv, "w") do io
        header = "name,nBin,nInt,nCont," *
                  join(["rootTime_run$k" for k in 1:n], ",") *
                  ",mean,min,max,nTimelimit"
        println(io, header)

        for name in allNames
            times = Float64[]
            statuses = String[]
            nBin = nInt = nCont = "?"
            for r in reports
                if haskey(r, name)
                    row = r[name]
                    push!(times, row.rootTime)
                    push!(statuses, row.status)
                    nBin, nInt, nCont = row.nBin, row.nInt, row.nCont
                else
                    push!(times, NaN)
                    push!(statuses, "missing")
                end
            end

            validTimes = filter(!isnan, times)
            meanT = isempty(validTimes) ? NaN : sum(validTimes) / length(validTimes)
            minT = isempty(validTimes) ? NaN : minimum(validTimes)
            maxT = isempty(validTimes) ? NaN : maximum(validTimes)
            nTimelimit = count(==("timelimit"), statuses)

            timesStr = join([isnan(t) ? "NA" : string(round(t, digits=3)) for t in times], ",")
            println(io, "$name,$nBin,$nInt,$nCont,$timesStr,$(round(meanT,digits=3)),$(round(minT,digits=3)),$(round(maxT,digits=3)),$nTimelimit")
        end
    end

    println("Aggregated $n run reports across $(length(allNames)) instances -> $outCsv")
end

main()
