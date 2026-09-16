#!/bin/bash

# Build a per-instance x per-config CSV: for every instance, which of the given configs
# solved it, and the heuristic time it took if so. Complements analyze_configs.sh (which
# aggregates per-config) with the instance-level view -- e.g. "does run4 solve instances
# baseline doesn't?"
#
# Usage:
#   ./instance_matrix.sh <folder1> [folder2] ...
#
# Example:
#   ./instance_matrix.sh fpfw_baseline fpfw_run1 fpfw_run2 fpfw_run3 fpfw_run4 fpfw_run5

if [ $# -lt 1 ]; then
    echo "Usage: ./instance_matrix.sh <folder1> [folder2] ..." >&2
    exit 1
fi

PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
COMP_RESULT="$PROJECT_DIR/compResult"

if ! command -v jq > /dev/null 2>&1; then
    echo "Error: jq not found; required to parse results.json" >&2
    exit 1
fi

for folder in "$@"; do
    if [ ! -d "$COMP_RESULT/$folder" ]; then
        echo "Error: $COMP_RESULT/$folder not found" >&2
        exit 1
    fi
done

# Warn (don't fail) if the given folders were run against different instance sets -
# e.g. miplib_selected was changed between runs. Folders from before
# instance_manifest.txt existed have no manifest and can't be verified.
FIRST_MANIFEST="$COMP_RESULT/$1/instance_manifest.txt"
if [ -f "$FIRST_MANIFEST" ]; then
    for folder in "$@"; do
        m="$COMP_RESULT/$folder/instance_manifest.txt"
        if [ ! -f "$m" ]; then
            echo "Warning: $folder has no recorded instance_manifest.txt (older run) - cannot verify it used the same instance set as $1" >&2
        elif ! diff -q "$FIRST_MANIFEST" "$m" > /dev/null 2>&1; then
            echo "Warning: $folder was run against a DIFFERENT instance set than $1 (instance_manifest.txt differs) - matrix rows below may not line up" >&2
        fi
    done
fi

OUT_DIR="$COMP_RESULT/config_comparison"
mkdir -p "$OUT_DIR"
CSV="$OUT_DIR/instance_matrix.csv"

# Instance list comes from the first folder given; every folder is expected to share the
# same fixed instance set (miplib_selected), since that's what submit_experiment.sh uses.
FIRST_FOLDER="$1"
mapfile -t instances < <(
    for d in "$COMP_RESULT/$FIRST_FOLDER"/*/; do
        basename "${d%/}"
    done | sort
)

# Header
{
    printf 'instance'
    for folder in "$@"; do
        printf ',%s_solved,%s_heurTime' "$folder" "$folder"
    done
    printf '\n'
} > "$CSV"

for instance in "${instances[@]}"; do
    row="$instance"
    for folder in "$@"; do
        results_json="$COMP_RESULT/$folder/$instance/results.json"
        if [ -f "$results_json" ]; then
            solved=$(jq -r '.solutionFound' "$results_json")
            if [ "$solved" = "true" ]; then
                heur_time=$(awk -v t="$(jq -r '.heurTime' "$results_json")" 'BEGIN { printf "%.2f", t }')
                row="${row},TRUE,${heur_time}"
            else
                row="${row},FALSE,"
            fi
        else
            row="${row},," # no results.json at all -- crashed/never ran
        fi
    done
    echo "$row" >> "$CSV"
done

echo "Instance matrix saved to $CSV"
