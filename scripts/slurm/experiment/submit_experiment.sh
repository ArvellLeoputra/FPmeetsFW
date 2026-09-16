#!/bin/bash

# Submit an FP-FW experiment run over every instance in INSTANCE_DIR, using a given
# settings/<cfgName>.cfg, into compResult/<folderName>/<instance>/{slurm_job.out,slurm_job.err,results.json}
#
# Usage:
#   ./submit_experiment.sh <cfgName> <folderName> [slurmWalltime] [seed] [throttle]
#
#   slurmWalltime  value for "#SBATCH --time" (default 1:00:00); NOT the heuristic
#                  time limit, which stays in the .cfg (timeLimit=...).
#   seed           optional; when given, the archived config's "seed=" line is
#                  overridden with this value and the result folder is suffixed
#                  "_s<seed>" so multiple seeds of the same config live side by
#                  side (compResult/<folderName>_s<seed>/). Also settable via the
#                  SEED env var.
#   throttle       max concurrent array tasks, i.e. the "%N" in --array=1-N%throttle
#                  (default 8). Raise it on a quiet cluster; same knob as -j in
#                  submit_sweep.sh.
#
# Example:
#   ./submit_experiment.sh fpfw fpfw_baseline
# Example:
#   ./submit_experiment.sh run3 fpfw_run3
# Example (5 seeds of run4 -> fpfw_run4_s1 .. fpfw_run4_s5):
#   for s in 1 2 3 4 5; do ./submit_experiment.sh run4 fpfw_run4 1:00:00 "$s"; done
# Example (32-wide instead of the default 8):
#   ./submit_experiment.sh run3 fpfw_run3 1:00:00 "" 32

if [ $# -lt 2 ]; then
    echo "Usage: ./submit_experiment.sh <cfgName> <folderName> [slurmWalltime] [seed] [throttle]" >&2
    exit 1
fi

# Positional args
CFG_NAME="$1"
FOLDER="$2"
TIME_LIMIT="${3:-1:00:00}"
SEED="${4:-${SEED:-}}"
THROTTLE="${5:-8}"

case "$THROTTLE" in ''|*[!0-9]*) echo "Error: throttle must be a positive integer, got '$THROTTLE'" >&2; exit 1 ;; esac
[ "$THROTTLE" -ge 1 ] || { echo "Error: throttle must be >= 1" >&2; exit 1; }

# When a seed is requested, give this run its own result tree so seeds don't
# overwrite each other. analyze_configs.sh globs "<folderName>_s*" to aggregate.
if [ -n "$SEED" ]; then
    case "$SEED" in
        ''|*[!0-9-]*) echo "Error: seed must be an integer, got: $SEED" >&2; exit 1 ;;
    esac
    FOLDER="${FOLDER}_s${SEED}"
fi

# Machine-specific paths, overridable via env (e.g. PROJECT_DIR=... ./submit_experiment.sh ...)
PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/instances/miplib_selected"
COMP_RESULT="$PROJECT_DIR/compResult"
CONFIG="$FPFW_DIR/settings/${CFG_NAME}.cfg"

# Check existence of config file
if [ ! -f "$CONFIG" ]; then
    echo "Error: Config not found: $CONFIG" >&2
    exit 1
fi

# Count instances up front; this becomes the SLURM array size
NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR" >&2
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

# (Re)create a clean result tree for this run. Each instance gets its own
# subdirectory holding slurm_job.out, slurm_job.err, and results.json
RESULT_DIR="$COMP_RESULT/$FOLDER"
rm -rf "$RESULT_DIR"
mkdir -p "$RESULT_DIR"

# Copy config here so results stay reproducible
RUN_CONFIG="$RESULT_DIR/config.cfg"
cp "$CONFIG" "$RUN_CONFIG"

# Override the seed in the archived config (not the source settings file) so
# every seed of a config is reproducible from its own compResult folder.
if [ -n "$SEED" ]; then
    if grep -q '^seed=' "$RUN_CONFIG"; then
        sed -i "s/^seed=.*/seed=$SEED/" "$RUN_CONFIG"
    else
        echo "seed=$SEED" >> "$RUN_CONFIG"
    fi
    echo "Seed override: seed=$SEED (folder $FOLDER)"
fi

# Write each instance's id/path once, so task IDs stay stable even if INSTANCE_DIR changes later
EXPERIMENT_LIST="$RESULT_DIR/experiment_list.tsv"
: > "$EXPERIMENT_LIST"
i=1
while IFS= read -r instance; do
    printf '%d\t%s\t%s\n' "$i" "$instance" "$INSTANCE_DIR/$instance" >> "$EXPERIMENT_LIST"
    i=$((i + 1))
done < <(ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort)

# Record the exact instance set this run used (sorted basenames), so later
# cross-run comparisons (analyze_configs.sh / instance_matrix.sh) can detect if
# INSTANCE_DIR's contents changed between runs instead of silently comparing
# mismatched instance sets.
cut -f2 "$EXPERIMENT_LIST" | sort > "$RESULT_DIR/instance_manifest.txt"

# Check envsubst is installed
if ! command -v envsubst > /dev/null 2>&1; then
    echo "Error: envsubst not found (part of gettext); required to render job_template.sh" >&2
    exit 1
fi

# Fill in the job template with submission-time values only
TEMPLATE="$(dirname "$0")/job_template.sh"
JOB_SCRIPT="$RESULT_DIR/job_script.sh"

export FOLDER TIME_LIMIT NUM_INSTANCES EXPERIMENT_LIST RESULT_DIR RUN_CONFIG FPFW_DIR THROTTLE
envsubst '$FOLDER,$TIME_LIMIT,$NUM_INSTANCES,$EXPERIMENT_LIST,$RESULT_DIR,$RUN_CONFIG,$FPFW_DIR,$THROTTLE' \
    < "$TEMPLATE" > "$JOB_SCRIPT"

# Submit the saved job script
sbatch "$JOB_SCRIPT" || { echo "ERROR: sbatch failed for $FOLDER"; exit 1; }

echo "Submitted: $FOLDER (${NUM_INSTANCES} instances, cfg=${CFG_NAME}${SEED:+, seed=${SEED}}, walltime=${TIME_LIMIT}, throttle=${THROTTLE})"
