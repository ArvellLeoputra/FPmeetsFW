using Plots

function plotIterates(
    scip::Ptr{SCIP.SCIP_},
    lpRows::Vector{Ptr{SCIP.SCIP_ROW}},
    lpCols::Vector{Ptr{SCIP.SCIP_COL}},
    colDict::Dict{Ptr{SCIP.SCIP_COL}, Int},
    intIdx::Vector{Int},
    xIterates::Vector{Vector{Float64}},
    xRoundIterates::Vector{Vector{Float64}},
    config::FPFWConfig,
    instanceName::String
)
    @assert length(intIdx) == 2 "plotIterates only supports 2D problems"

    i1, i2 = intIdx[1], intIdx[2]

    p = plot(
        xlims=(-0.1, 1.1), ylims=(-0.1, 1.1),
        xlabel="x[$(i1)]", ylabel="x[$(i2)]",
        title="$instanceName | $(config.norm) | $(config.fwVariant) | $(config.lineSearch)",
        legend=:topright,
        aspect_ratio=:equal
    )

    # Draw [0,1]² boundary
    plot!(p, [0,1,1,0,0], [0,0,1,1,0], color=:black, linewidth=1, label=false)

    # Draw constraint lines for rows that only involve the two integer variables
    t_range = LinRange(-0.1, 1.1, 200)
    for row in lpRows
        nnz = SCIP.SCIProwGetNNonz(row)
        cols_row = unsafe_wrap(Array, SCIP.SCIProwGetCols(row), nnz)
        vals_row = unsafe_wrap(Array, SCIP.SCIProwGetVals(row), nnz)

        a1, a2 = 0.0, 0.0
        pure2D = true
        for k in 1:nnz
            idx = get(colDict, cols_row[k], 0)
            if idx == i1
                a1 = vals_row[k]
            elseif idx == i2
                a2 = vals_row[k]
            else
                pure2D = false
                break
            end
        end
        !pure2D && continue

        lhs = SCIP.SCIProwGetLhs(row)
        rhs = SCIP.SCIProwGetRhs(row)

        for b in (lhs, rhs)
            abs(b) > 1e29 && continue
            if abs(a2) > 1e-10
                y_vals = [(b - a1 * t) / a2 for t in t_range]
                plot!(p, collect(t_range), y_vals, linestyle=:dot, color=:gray, linewidth=1, label=false)
            elseif abs(a1) > 1e-10
                vline!(p, [b / a1], linestyle=:dot, color=:gray, linewidth=1, label=false)
            end
        end
    end

    # Mark integer feasible points in {0,1}²
    feasiblePts = Tuple{Float64,Float64}[]
    for v1 in [0.0, 1.0], v2 in [0.0, 1.0]
        sol = zeros(Float64, length(lpCols))
        sol[i1] = v1
        sol[i2] = v2
        if isSolutionLPFeasible(scip, lpRows, lpCols, sol, colDict)
            push!(feasiblePts, (v1, v2))
        end
    end
    if !isempty(feasiblePts)
        scatter!(p, [pt[1] for pt in feasiblePts], [pt[2] for pt in feasiblePts],
            marker=:star5, markersize=12, color=:green, label="Integer feasible")
    end

    # Plot x iterates trajectory
    xs = [x[i1] for x in xIterates]
    ys = [x[i2] for x in xIterates]
    plot!(p, xs, ys, color=:blue, linewidth=1, label=false)
    scatter!(p, xs, ys, marker=:circle, markersize=6, color=:blue, label="x")
    for (k, (x1, x2)) in enumerate(zip(xs, ys))
        annotate!(p, x1 + 0.02, x2 + 0.02, text("$k", 8, :blue))
    end

    # Plot xRound iterates
    xrs = [xr[i1] for xr in xRoundIterates]
    yrs = [xr[i2] for xr in xRoundIterates]
    scatter!(p, xrs, yrs, marker=:rect, markersize=6, color=:red, label="xRound")

    # Save to plots/ folder
    mkpath("plots")
    filename = "plots/$(instanceName)_$(config.norm)_$(config.fwVariant)_$(config.lineSearch).png"
    savefig(p, filename)
    println("Plot saved to $filename")
end