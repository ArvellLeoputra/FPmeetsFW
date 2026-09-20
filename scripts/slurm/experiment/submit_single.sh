#!/bin/bash
# Runs exactly one instance as a single SLURM job, non-exclusive and writes results.json into compResult/<folder>/<instance>/
# Created for backfilling instances an array job never finished (e.g. silently killed by the wall-clock limit)
#
# Usage:
#   ./submit_single.sh [-R] <instance_basename> <folder> [slurmWalltime]
#
#   -R             run via clean.jl instead of main.jl (a discarded warm-up
#                  solve, then the real one). Match this to whether the
#                  original sweep used -R, or timing won't be comparable.
#   slurmWalltime  value for "#SBATCH --time" (default 1:00:00). Worth raising
#                  above the array job's original walltime when re-running an
#                  instance that was killed for running out of time.
#
# Example (re-run rmatr200-p5 into the already-existing fpfw_run2 folder):
#   ./submit_single.sh rmatr200-p5 fpfw_run2
# Example (same, but with a longer walltime since it was killed by the original one):
#   ./submit_single.sh rmatr200-p5 fpfw_run2 4:00:00
# Example (original sweep used -R, so this backfill must too):
#   ./submit_single.sh -R rmatr200-p5 fpfw_run2

set -euo pipefail

USE_CLEAN=0
while getopts "R" opt; do
    case "$opt" in
        R) USE_CLEAN=1 ;;
        *) echo "Usage: ./submit_single.sh [-R] <instance_basename> <folder> [slurmWalltime]" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -lt 2 ]; then
    echo "Usage: ./submit_single.sh [-R] <instance_basename> <folder> [slurmWalltime]" >&2
    exit 1
fi

INSTANCE_BASENAME="$1"
FOLDER="$2"
TIME_LIMIT="${3:-1:00:00}"

RUNNER="main.jl"
[ "$USE_CLEAN" = 1 ] && RUNNER="clean.jl"

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
#SBATCH --time=$TIME_LIMIT
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --output=$INSTANCE_RESULT_DIR/slurm_job.out
#SBATCH --error=$INSTANCE_RESULT_DIR/slurm_job.err

export PATH="/home/htc/aleoputra/scratch/julia-1.12.6/bin:\$PATH"

echo "Running instance: $INSTANCE_PATH"
echo "Config: $CONFIG"
echo "Node: \$(hostname)"

julia --project=$FPFW_DIR $FPFW_DIR/$RUNNER "$INSTANCE_PATH" "$CONFIG" "resultsDir=$INSTANCE_RESULT_DIR"
EOF

echo "Submitted single run: $INSTANCE_BASENAME -> $FOLDER (walltime=$TIME_LIMIT, runner=$RUNNER)"
