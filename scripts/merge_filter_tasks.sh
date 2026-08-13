#!/bin/bash
# Merges the per-instance CSV rows written by SLURM array tasks
# (scripts/filter_instances.jl ... task <taskIndex>) into a single report,
# with one shared header line.
#
# Usage:
#   bash scripts/merge_filter_tasks.sh <taskDir> <outCsv>
#
# Example:
#   bash scripts/merge_filter_tasks.sh misc/filter_tasks misc/filter_report.csv

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: bash scripts/merge_filter_tasks.sh <taskDir> <outCsv>" >&2
    exit 1
fi

taskDir="$1"
outCsv="$2"

taskFiles=("$taskDir"/task_*.csv)
if [ ! -e "${taskFiles[0]}" ]; then
    echo "No task_*.csv files found in $taskDir" >&2
    exit 1
fi

# Header from the first task file, then every task file's single data row (line 2).
head -n 1 "${taskFiles[0]}" > "$outCsv"
for f in "$taskDir"/task_*.csv; do
    sed -n '2p' "$f" >> "$outCsv"
done

nRows=$(($(wc -l < "$outCsv") - 1))
echo "Merged $nRows task rows into $outCsv"
