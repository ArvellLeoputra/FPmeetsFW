#!/bin/bash
# SLURM array job: one task per MIPLIB instance, running the same root-node-only
# check as `julia --project scripts/filter_instances.jl <dir> ... task <taskIndex>`
# (presolve on, node_limit=1, no heuristics/cuts)
#
# Usage (from any directory):
#   sbatch /home/htc/aleoputra/project/FPmeetsFW/scripts/slurm/instance_filter/filter_array.sh
#
# After all tasks finish, merge the per-task CSVs into one report:
#   bash scripts/merge_filter.sh misc/filter_tasks misc/filter_report.csv

#SBATCH --job-name=fpfw-filter
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --array=1-240
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/htc/aleoputra/project/FPmeetsFW/misc/slurm_logs/filter_%A_%a.out
#SBATCH --error=/home/htc/aleoputra/project/FPmeetsFW/misc/slurm_logs/filter_%A_%a.err

set -euo pipefail

export PATH="/home/htc/aleoputra/scratch/julia-1.12.6/bin:$PATH"

FPFW_DIR="/home/htc/aleoputra/project/FPmeetsFW"
cd "$FPFW_DIR"
mkdir -p misc/slurm_logs misc/filter_tasks

INSTANCE_DIR="/home/htc/aleoputra/project/instances/miplib_full"
ROOT_TIME_LIMIT=300.0
OUT_DIR="misc/filter_tasks"

julia --project scripts/filter_instances.jl "$INSTANCE_DIR" "$OUT_DIR" "$ROOT_TIME_LIMIT" task "$SLURM_ARRAY_TASK_ID"