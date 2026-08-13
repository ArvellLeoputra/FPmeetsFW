#!/bin/bash

# Wipe filter_array.sh's output before a fresh submission (stale logs/CSVs from
# a previous/failed array run). Run this ONCE by hand before `sbatch filter_array.sh`
# — never from inside the array job itself, since with 240 concurrent tasks a
# per-task rm -rf would race and delete other tasks' in-progress/finished output.
#
# Usage: ./reset_filter.sh

FPFW_DIR="/home/htc/aleoputra/project/FPmeetsFW"

rm -rf "$FPFW_DIR/misc/slurm_logs" "$FPFW_DIR/misc/filter_tasks"
mkdir -p "$FPFW_DIR/misc/slurm_logs" "$FPFW_DIR/misc/filter_tasks"

echo "Reset $FPFW_DIR/misc/slurm_logs and $FPFW_DIR/misc/filter_tasks"
