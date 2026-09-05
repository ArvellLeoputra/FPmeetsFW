#!/bin/bash

# Compare aggregate stats across multiple FP-FW config runs, side by side.
#
# Each argument is a config name. If per-seed result folders "<config>_s<seed>"
# exist under compResult/ (as produced by submit_experiment.sh ... <seed>), the
# script aggregates them: it reports mean / sample SD / min / max of the solved
# count across seeds, plus per-seed means of the other columns. If no seed
# folders are found it falls back to a single folder named exactly <config>
# (the previous behaviour), reported as SEEDS=1, SD=0.
#
# For each folder it reads compResult/config_comparison/<folder>/solution_summary.txt,
# generating it via analyze_results.sh first if missing.
#
# Writes into compResult/config_comparison/:
#   config_comparison.txt          leaderboard (one row per config)
#   config_comparison.csv          one row per config, with SD/min/max
#   config_comparison_by_seed.csv  one row per (config, seed) - use this for
#                                  per-seed significance tests / plots
#
# Usage:
#   ./analyze_configs.sh <config1> [config2] ...
#
# Examples:
#   ./analyze_configs.sh fpfw_baseline fpfw_run1 fpfw_run2 fpfw_run3 fpfw_run4 fpfw_run5
#   # with fpfw_run4_s1 .. fpfw_run4_s5 present, this aggregates them:
#   ./analyze_configs.sh fpfw_run4

if [ $# -lt 1 ]; then
    echo "Usage: ./analyze_configs.sh <config1> [config2] ..." >&2
    exit 1
fi

PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
COMP_RESULT="$PROJECT_DIR/compResult"
SCRIPT_DIR="$(dirname "$0")"

if ! command -v jq > /dev/null 2>&1; then
    echo "Error: jq not found; required (via analyze_results.sh) to parse results.json" >&2
    exit 1
fi

OUT_DIR="$COMP_RESULT/config_comparison"
mkdir -p "$OUT_DIR"
TABLE="$OUT_DIR/config_comparison.txt"
CSV="$OUT_DIR/config_comparison.csv"
BYSEED_CSV="$OUT_DIR/config_comparison_by_seed.csv"

# Extract one labeled numeric field from a solution_summary.txt, e.g.
#   field summary.txt "Total instances:"   -> "208"
field() {
    local file="$1" label="$2"
    grep -F "$label" "$file" | head -1 | sed -E "s/^[[:space:]]*${label//\//\\/}[[:space:]]*//" | awk '{print $1}'
}

# Ensure compResult/config_comparison/<folder>/solution_summary.txt exists.
# Returns 0 on success. analyze_results.sh reads the folder's archived
# config.cfg, so the (unused-here) cfgName arg is just the folder name.
ensure_summary() {
    local folder="$1"
    local summary="$COMP_RESULT/config_comparison/$folder/solution_summary.txt"
    if [ ! -f "$summary" ]; then
        echo "No summary yet for $folder, generating via analyze_results.sh..." >&2
        ( cd "$SCRIPT_DIR" && ./analyze_results.sh "$folder" "$folder" > /dev/null )
    fi
    [ -f "$summary" ]
}

# mean / sample-SD / min / max of numbers read from stdin (one per line).
# Prints "mean sd min max"; SD is 0 when n < 2.
statline() {
    awk '
        { v = $1 + 0; x[NR] = v; s += v
          if (NR == 1 || v < mn) mn = v
          if (NR == 1 || v > mx) mx = v }
        END {
            if (NR == 0) { print "0 0 0 0"; exit }
            m = s / NR
            for (i = 1; i <= NR; i++) d += (x[i] - m) * (x[i] - m)
            sd = (NR > 1) ? sqrt(d / (NR - 1)) : 0
            printf "%.2f %.2f %d %d\n", m, sd, mn, mx
        }'
}

# mean of numbers read from stdin (one per line).
meanline() {
    awk '{ s += $1; n++ } END { if (n == 0) print 0; else printf "%.2f\n", s / n }'
}

echo "config,seed,folder,solved,solve_pct,failed,time_limit,scip_time_limit,iter_limit,fw_infeasible,unknown,total,avg_total_time_s,avg_heur_time_s" > "$BYSEED_CSV"

declare -a rows=()

for config in "$@"; do
    # Decide grouped (per-seed) vs single folder.
    seed_folders=()
    for d in "$COMP_RESULT/${config}"_s*/; do
        [ -d "$d" ] || continue
        seed_folders+=("$(basename "${d%/}")")
    done

    if [ ${#seed_folders[@]} -gt 0 ]; then
        folders=("${seed_folders[@]}")
    elif [ -d "$COMP_RESULT/$config" ]; then
        folders=("$config")
    else
        echo "Warning: no folders for '$config' (looked for $config and ${config}_s*), skipping" >&2
        continue
    fi

    solved_list=""; failed_list=""; tl_list=""; scip_list=""; iter_list=""
    fwinf_list=""; unk_list=""; time_list=""; heur_list=""; pct_list=""
    total_val=""; nseeds=0

    for f in "${folders[@]}"; do
        if ! ensure_summary "$f"; then
            echo "Warning: could not generate summary for $f, skipping" >&2
            continue
        fi
        s="$COMP_RESULT/config_comparison/$f/solution_summary.txt"

        seed="${f##*_s}"
        [ "$seed" = "$f" ] && seed="-"

        sv=$(field "$s" "Solutions found:")
        pc=$(grep -F "Solutions found:" "$s" | head -1 | grep -oE '\([0-9.]+%\)' | tr -d '(%)')
        fl=$(field "$s" "Failed:")
        tl=$(field "$s" "  Time limit:")
        stl=$(field "$s" "  SCIP time limit:")
        il=$(field "$s" "  Iteration limit:")
        fi_=$(field "$s" "  FW infeasible:")
        un=$(field "$s" "  Unknown:")
        at=$(field "$s" "Average total time:" | sed 's/s$//')
        ah=$(field "$s" "Average heur time:" | sed 's/s$//')
        tot=$(field "$s" "Total instances:")
        [ -n "$tot" ] && total_val="$tot"

        nseeds=$((nseeds + 1))
        solved_list+="${sv:-0}"$'\n'
        pct_list+="${pc:-0}"$'\n'
        failed_list+="${fl:-0}"$'\n'
        tl_list+="${tl:-0}"$'\n'
        scip_list+="${stl:-0}"$'\n'
        iter_list+="${il:-0}"$'\n'
        fwinf_list+="${fi_:-0}"$'\n'
        unk_list+="${un:-0}"$'\n'
        time_list+="${at:-0}"$'\n'
        heur_list+="${ah:-0}"$'\n'

        echo "$config,$seed,$f,${sv:-0},${pc:-0},${fl:-0},${tl:-0},${stl:-0},${il:-0},${fi_:-0},${un:-0},${tot:-0},${at:-0},${ah:-0}" >> "$BYSEED_CSV"
    done

    [ "$nseeds" -eq 0 ] && continue

    read -r sv_mean sv_sd sv_min sv_max <<< "$(printf '%s' "$solved_list" | statline)"
    pc_mean=$(printf '%s' "$pct_list"    | meanline)
    fl_mean=$(printf '%s' "$failed_list" | meanline)
    tl_mean=$(printf '%s' "$tl_list"     | meanline)
    stl_mean=$(printf '%s' "$scip_list"  | meanline)
    il_mean=$(printf '%s' "$iter_list"   | meanline)
    fi_mean=$(printf '%s' "$fwinf_list"  | meanline)
    un_mean=$(printf '%s' "$unk_list"    | meanline)
    at_mean=$(printf '%s' "$time_list"   | meanline)
    ah_mean=$(printf '%s' "$heur_list"   | meanline)

    # sort key | config | seeds | solved_mean | solved_sd | solved_min | solved_max
    #   | solve%_mean | failed_mean | tl_mean | scip_tl_mean | iter_mean
    #   | fw_infeas_mean | unknown_mean | total | avg_time_mean | avg_heur_mean
    rows+=("${sv_mean}|${config}|${nseeds}|${sv_mean}|${sv_sd}|${sv_min}|${sv_max}|${pc_mean}|${fl_mean}|${tl_mean}|${stl_mean}|${il_mean}|${fi_mean}|${un_mean}|${total_val:-0}|${at_mean}|${ah_mean}")
done

if [ ${#rows[@]} -eq 0 ]; then
    echo "Error: no valid config summaries found" >&2
    exit 1
fi

# Best mean solve count first
IFS=$'\n' sorted=($(printf '%s\n' "${rows[@]}" | sort -t'|' -k1,1rn))
unset IFS

{
echo "Config comparison (per-seed aggregate) - $(date)"
echo "============================================================================================================================================"
printf "%-16s %-6s %-9s %-7s %-6s %-6s %-8s %-9s %-7s %-7s %-8s %-7s %-9s %-7s %-8s %-8s\n" \
    "CONFIG" "SEEDS" "SOLVED~" "SD" "MIN" "MAX" "SOLVE%~" "FAILED~" "TL~" "SCIPTL~" "ITERLIM~" "FWINF~" "UNKNOWN~" "TOTAL" "AVGTIME~" "AVGHEUR~"
echo "--------------------------------------------------------------------------------------------------------------------------------------------"
for row in "${sorted[@]}"; do
    IFS='|' read -r _ config seeds sv_mean sv_sd sv_min sv_max pc_mean fl_mean tl_mean stl_mean il_mean fi_mean un_mean total at_mean ah_mean <<< "$row"
    printf "%-16s %-6s %-9s %-7s %-6s %-6s %-8s %-9s %-7s %-7s %-8s %-7s %-9s %-7s %-8s %-8s\n" \
        "$config" "$seeds" "$sv_mean" "$sv_sd" "$sv_min" "$sv_max" "${pc_mean}%" "$fl_mean" "$tl_mean" "$stl_mean" "$il_mean" "$fi_mean" "$un_mean" "$total" "${at_mean}s" "${ah_mean}s"
done
echo ""
echo "~ = mean across seeds. SD is the sample standard deviation of SOLVED across seeds (0 when SEEDS=1)."
echo "Per-(config,seed) rows: $BYSEED_CSV"
} | tee "$TABLE"

{
echo "config,seeds,solved_mean,solved_sd,solved_min,solved_max,solve_pct_mean,failed_mean,time_limit_mean,scip_time_limit_mean,iter_limit_mean,fw_infeasible_mean,unknown_mean,total,avg_total_time_s_mean,avg_heur_time_s_mean"
for row in "${sorted[@]}"; do
    IFS='|' read -r _ config seeds sv_mean sv_sd sv_min sv_max pc_mean fl_mean tl_mean stl_mean il_mean fi_mean un_mean total at_mean ah_mean <<< "$row"
    echo "$config,$seeds,$sv_mean,$sv_sd,$sv_min,$sv_max,$pc_mean,$fl_mean,$tl_mean,$stl_mean,$il_mean,$fi_mean,$un_mean,$total,$at_mean,$ah_mean"
done
} > "$CSV"

echo ""
echo "Comparison saved to:"
echo "  $TABLE"
echo "  $CSV"
echo "  $BYSEED_CSV"
