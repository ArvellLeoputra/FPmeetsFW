#!/bin/bash
#
# submit_packed.sh - run one config over every instance in miplib_selected on a
# SINGLE, exclusively-reserved node, packing up to <concurrency> instances at a
# time - instead of relying on SLURM's array scheduler (submit_experiment.sh),
# which can place tasks on any node shared with other users' jobs.
#
# This reserves exactly one whole node for the entire run (--exclusive) - no
# other job can land on it while it's running - while still packing many of
# YOUR OWN instances onto it concurrently, so you're not wasting the other
# ~100+ idle cores like a plain --exclusive array task would (see
# submit_sweep.sh's -x, which reserves one node PER TASK, not one node total).
#
# Tradeoff vs submit_experiment.sh: bounded by one node's capacity (RAM-bound
# around ~60 instances at 16G each on a ~1TB node - stay well under that with
# <concurrency>), and instances now run in sequential waves of <concurrency>
# rather than being spread across the whole cluster, so total wall-clock is
# roughly ceil(numInstances/concurrency) x typical-instance-time. In exchange,
# every instance in this run shares hardware with nothing but itself and its
# own siblings - no cross-user memory-bandwidth/cache/turbo noise.
#
# Usage:
#   ./submit_packed.sh <cfgName> <folderName> [concurrency] [slurmWalltime]
#
#   concurrency    max instances running at once on the one reserved node
#                  (default 30).
#   slurmWalltime  value for "#SBATCH --time" for the WHOLE run (default
#                  4:00:00) - covers ceil(numInstances/concurrency) waves,
#                  not just one instance, since it's one job holding the node
#                  for the entire run.
#
# Results land exactly like submit_experiment.sh's:
#   compResult/<folderName>/<instance>/{slurm_job.out,slurm_job.err,results.json}
# so analyze_results.sh / analyze_configs.sh / instance_matrix.sh work unchanged.
#
# Example:
#   ./submit_packed.sh fpfw fpfw_packed 30
# Example (longer walltime, lower concurrency for a heavier config):
#   ./submit_packed.sh run3 fpfw_run3_packed 20 8:00:00

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: ./submit_packed.sh <cfgName> <folderName> [concurrency] [slurmWalltime]" >&2
    exit 1
fi

CFG_NAME="$1"
FOLDER="$2"
CONCURRENCY="${3:-30}"
TIME_LIMIT="${4:-4:00:00}"

case "$CONCURRENCY" in ''|*[!0-9]*) echo "Error: concurrency must be a positive integer, got '$CONCURRENCY'" >&2; exit 1 ;; esac
[ "$CONCURRENCY" -ge 1 ] || { echo "Error: concurrency must be >= 1" >&2; exit 1; }

# Machine-specific paths, overridable via env (same convention as the other submit_*.sh)
PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/instances/miplib_selected"
COMP_RESULT="$PROJECT_DIR/compResult"
CONFIG="$FPFW_DIR/settings/${CFG_NAME}.cfg"
JULIA_BIN="${JULIA_BIN:-/home/htc/aleoputra/scratch/julia-1.12.6/bin}"

if [ ! -f "$CONFIG" ]; then
    echo "Error: Config not found: $CONFIG" >&2
    exit 1
fi

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR" >&2
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

# (Re)create a clean result tree for this run.
RESULT_DIR="$COMP_RESULT/$FOLDER"
rm -rf "$RESULT_DIR"
mkdir -p "$RESULT_DIR"

# Copy config here so results stay reproducible
RUN_CONFIG="$RESULT_DIR/config.cfg"
cp "$CONFIG" "$RUN_CONFIG"

# Same manifest convention as submit_experiment.sh/submit_sweep.sh - full
# filenames, sorted - so analyze_configs.sh / instance_matrix.sh can detect if
# this run used a different instance set than one submitted the other way.
# Doubles as the work list handed to xargs below.
INSTANCE_LIST="$RESULT_DIR/instance_manifest.txt"
ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort > "$INSTANCE_LIST"

# Render the one-node packed job script.
JOB_SCRIPT="$RESULT_DIR/job_script.sh"
cat > "$JOB_SCRIPT" <<TEMPLATE
#!/bin/bash
#SBATCH --job-name=$FOLDER
#SBATCH --time=$TIME_LIMIT
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --output=$RESULT_DIR/slurm_job.out
#SBATCH --error=$RESULT_DIR/slurm_job.err

set -euo pipefail
export PATH="$JULIA_BIN:\$PATH"

echo "Packed run on \$(hostname): $NUM_INSTANCES instances, concurrency=$CONCURRENCY"
echo "Config: $RUN_CONFIG"

run_one() {
    local instance="\$1"
    local base="\${instance%.mps.gz}"
    base="\${base%.mps}"
    local instance_result_dir="$RESULT_DIR/\$base"
    mkdir -p "\$instance_result_dir"
    julia --project=$FPFW_DIR $FPFW_DIR/main.jl "$INSTANCE_DIR/\$instance" "$RUN_CONFIG" "resultsDir=\$instance_result_dir" \
        > "\$instance_result_dir/slurm_job.out" 2> "\$instance_result_dir/slurm_job.err"
}
export -f run_one

cat "$INSTANCE_LIST" | xargs -P $CONCURRENCY -I{} bash -c 'run_one "\$@"' _ {}
TEMPLATE
chmod +x "$JOB_SCRIPT"

sbatch "$JOB_SCRIPT" || { echo "ERROR: sbatch failed for $FOLDER"; exit 1; }

echo "Submitted: $FOLDER ($NUM_INSTANCES instances, cfg=${CFG_NAME}, concurrency=${CONCURRENCY}, 1 exclusive node, walltime=${TIME_LIMIT})"
