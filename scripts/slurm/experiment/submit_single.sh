#!/bin/bash
# Runs exactly one instance as a single SLURM job, non-exclusive and writes results.json into compResult/<folder>/<instance>/
# Created for backfilling instances an array job never finished (e.g. silently killed by the wall-clock limit)
#
# Usage:
#   ./submit_single.sh <instance_basename> <folder>
#
# Example (re-run rmatr200-p5 into the already-existing fpfw_run2 folder):
#   ./submit_single.sh rmatr200-p5 fpfw_run2

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: ./submit_single.sh <instance_basename> <folder>" >&2
    exit 1
fi

INSTANCE_BASENAME="$1"
FOLDER="$2"

PROJECT_DIR="/home/htc/aleoputra/project"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/instances/miplib_selected"
COMP_RESULT="$PROJECT_DIR/compResult"
RESULT_DIR="$COMP_RESULT/$FOLDER"
CONFIG="$RESULT_DIR/config.cfg"

INSTANCE_PATH="$INSTANCE_DIR/${INSTANCE_BASENAME}.mps.gz"
if [ ! -f "$INSTANCE_PATH" ]; then
    echo "Error: $INSTANCE_PATH not found" >&2
    exit 1
fi
if [ ! -f "$CONFIG" ]; then
    echo "Error: $CONFIG not found (expected the archived config from submit_experiment.sh)" >&2
    exit 1
fi

INSTANCE_RESULT_DIR="$RESULT_DIR/$INSTANCE_BASENAME"
mkdir -p "$INSTANCE_RESULT_DIR"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=single_${INSTANCE_BASENAME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --output=$INSTANCE_RESULT_DIR/slurm_job.out
#SBATCH --error=$INSTANCE_RESULT_DIR/slurm_job.err

export PATH="/home/htc/aleoputra/scratch/julia-1.12.6/bin:\$PATH"

echo "Running instance: $INSTANCE_PATH"
echo "Config: $CONFIG"
echo "Node: \$(hostname)"

julia --project=$FPFW_DIR $FPFW_DIR/main.jl "$INSTANCE_PATH" "$CONFIG" "resultsDir=$INSTANCE_RESULT_DIR"
EOF

echo "Submitted single run: $INSTANCE_BASENAME -> $FOLDER"
