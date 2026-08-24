#!/bin/bash
#SBATCH --job-name=$FOLDER
#SBATCH --time=$TIME_LIMIT
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --array=1-$NUM_INSTANCES%8
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

# Add Julia to PATH
export PATH="/home/htc/aleoputra/scratch/julia-1.12.6/bin:$PATH"

# Look up this task's instance
LINE=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" '$1 == id {print; exit}' "$EXPERIMENT_LIST")

if [ -z "$LINE" ]; then
    exit 1
fi

# Parse instance name and path
IFS=$'\t' read -r _ INSTANCE INSTANCE_PATH <<< "$LINE"

# Create a result directory for each instance
BASENAME=$(basename "$INSTANCE" .mps.gz)
BASENAME=$(basename "$BASENAME" .mps)
INSTANCE_RESULT_DIR="$RESULT_DIR/$BASENAME"
mkdir -p "$INSTANCE_RESULT_DIR"

# Redirect output to that instance's log files
exec > "$INSTANCE_RESULT_DIR/slurm_job.out" 2> "$INSTANCE_RESULT_DIR/slurm_job.err"

# Log run metadata
echo "Running instance: $INSTANCE_PATH"
echo "SLURM task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Config: $RUN_CONFIG"

# Run the FP-FW heuristic, writing results.json into this instance's folder.
julia --project=$FPFW_DIR $FPFW_DIR/main.jl "$INSTANCE_PATH" "$RUN_CONFIG" "resultsDir=$INSTANCE_RESULT_DIR"
