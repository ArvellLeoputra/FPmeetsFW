#!/bin/bash

# Compare aggregate stats across multiple FP-FW config runs, side by side.
# Reads each folder's compResult/<folder>/result/solution_summary.txt (generating
# it via analyze_results.sh first if missing), and writes one consolidated
# leaderboard into compResult/config_comparison/.
#
# Usage:
#   ./analyze_configs.sh <folder1> [folder2] [folder3] ...
#
# Example:
#   ./analyze_configs.sh fpfw_baseline fpfw_run1 fpfw_run2 fpfw_run3 fpfw_run4 fpfw_run5 fpfw_euc

if [ $# -lt 1 ]; then
    echo "Usage: ./analyze_configs.sh <folder1> [folder2] ..." >&2
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

# Extract one labeled numeric field from a solution_summary.txt, e.g.
#   field summary.txt "Total instances:"   -> "208"
field() {
    local file="$1" label="$2"
    grep -F "$label" "$file" | head -1 | sed -E "s/^[[:space:]]*${label//\//\\/}[[:space:]]*//" | awk '{print $1}'
}

declare -a rows=()

for folder in "$@"; do
    RESULT_DIR="$COMP_RESULT/$folder"
    SUMMARY="$COMP_RESULT/config_comparison/$folder/solution_summary.txt"

    if [ ! -d "$RESULT_DIR" ]; then
        echo "Warning: $RESULT_DIR not found, skipping" >&2
        continue
    fi

    if [ ! -f "$SUMMARY" ]; then
        echo "No summary yet for $folder, generating via analyze_results.sh..."
        ( cd "$SCRIPT_DIR" && ./analyze_results.sh "$folder" "$folder" > /dev/null )
    fi

    if [ ! -f "$SUMMARY" ]; then
        echo "Warning: could not generate summary for $folder, skipping" >&2
        continue
    fi

    total=$(field "$SUMMARY" "Total instances:")
    found=$(field "$SUMMARY" "Solutions found:")
    found_pct=$(grep -F "Solutions found:" "$SUMMARY" | head -1 | grep -oE '\([0-9.]+%\)' | tr -d '(%)')
    failed=$(field "$SUMMARY" "Failed:")
    tl=$(field "$SUMMARY" "  Time limit:")
    scip_tl=$(field "$SUMMARY" "  SCIP time limit:")
    iter_l=$(field "$SUMMARY" "  Iteration limit:")
    fw_infeas=$(field "$SUMMARY" "  FW infeasible:")
    unknown=$(field "$SUMMARY" "  Unknown:")
    avg_time=$(field "$SUMMARY" "Average total time:" | sed 's/s$//')
    avg_heur=$(field "$SUMMARY" "Average heur time:" | sed 's/s$//')

    rows+=("${found_pct:-0}|$folder|$total|$found|$failed|$tl|$scip_tl|$iter_l|$fw_infeas|$unknown|$avg_time|$avg_heur")
done

if [ ${#rows[@]} -eq 0 ]; then
    echo "Error: no valid config summaries found" >&2
    exit 1
fi

# Sort by solve rate, best first
IFS=$'\n' sorted=($(printf '%s\n' "${rows[@]}" | sort -t'|' -k1,1rn))
unset IFS

{
echo "Config comparison - $(date)"
echo "=========================================================================================================="
printf "%-16s %-6s %-8s %-8s %-6s %-6s %-8s %-6s %-8s %-8s %-8s %-8s\n" \
    "FOLDER" "SOLVED" "SOLVE%" "FAILED" "TL" "SCIPTL" "ITERLIM" "FWINF" "UNKNOWN" "TOTAL" "AVGTIME" "AVGHEUR"
echo "----------------------------------------------------------------------------------------------------------"
for row in "${sorted[@]}"; do
    IFS='|' read -r pct folder total found failed tl scip_tl iter_l fw_infeas unknown avg_time avg_heur <<< "$row"
    printf "%-16s %-6s %-8s %-8s %-6s %-6s %-8s %-6s %-8s %-8s %-8s %-8s\n" \
        "$folder" "$found" "${pct}%" "$failed" "$tl" "$scip_tl" "$iter_l" "$fw_infeas" "$unknown" "$total" "${avg_time}s" "${avg_heur}s"
done
} | tee "$TABLE"

{
echo "folder,solved,solve_pct,failed,time_limit,scip_time_limit,iter_limit,fw_infeasible,unknown,total,avg_total_time_s,avg_heur_time_s"
for row in "${sorted[@]}"; do
    IFS='|' read -r pct folder total found failed tl scip_tl iter_l fw_infeas unknown avg_time avg_heur <<< "$row"
    echo "$folder,$found,$pct,$failed,$tl,$scip_tl,$iter_l,$fw_infeas,$unknown,$total,$avg_time,$avg_heur"
done
} > "$CSV"

echo ""
echo "Comparison saved to $TABLE and $CSV"
