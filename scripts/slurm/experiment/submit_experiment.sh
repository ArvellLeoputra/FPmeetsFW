#!/bin/bash

# Submit an FP-FW experiment run over every instance in INSTANCE_DIR, using a given
# settings/<cfgName>.cfg, into compResult/<folderName>/<instance>/{slurm_job.out,slurm_job.err,results.json}
#
# Usage:
#   ./submit_experiment.sh <cfgName> <folderName> [timeLimit]
#
# Example (replaces the old submit_fpfw.sh):
#   ./submit_experiment.sh fpfw fpfw_baseline
# Example (replaces the old submit_run3.sh):
#   ./submit_experiment.sh run3 fpfw_run3

if [ $# -lt 2 ]; then
    echo "Usage: ./submit_experiment.sh <cfgName> <folderName> [timeLimit]" >&2
    exit 1
fi

# Positional args
CFG_NAME="$1"
FOLDER="$2"
TIME_LIMIT="${3:-30:00}"

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

# Write each instance's id/path once, so task IDs stay stable even if INSTANCE_DIR changes later
EXPERIMENT_LIST="$RESULT_DIR/experiment_list.tsv"
: > "$EXPERIMENT_LIST"
i=1
while IFS= read -r instance; do
    printf '%d\t%s\t%s\n' "$i" "$instance" "$INSTANCE_DIR/$instance" >> "$EXPERIMENT_LIST"
    i=$((i + 1))
done < <(ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort)

# Check envsubst is installed
if ! command -v envsubst > /dev/null 2>&1; then
    echo "Error: envsubst not found (part of gettext); required to render job_template.sh" >&2
    exit 1
fi

# Fill in the job template with submission-time values only
TEMPLATE="$(dirname "$0")/job_template.sh"
JOB_SCRIPT="$RESULT_DIR/job_script.sh"

export FOLDER TIME_LIMIT NUM_INSTANCES EXPERIMENT_LIST RESULT_DIR RUN_CONFIG FPFW_DIR
envsubst '$FOLDER,$TIME_LIMIT,$NUM_INSTANCES,$EXPERIMENT_LIST,$RESULT_DIR,$RUN_CONFIG,$FPFW_DIR' \
    < "$TEMPLATE" > "$JOB_SCRIPT"

# Submit the saved job script
sbatch "$JOB_SCRIPT" || { echo "ERROR: sbatch failed for $FOLDER"; exit 1; }

echo "Submitted: $FOLDER (${NUM_INSTANCES} instances, cfg=${CFG_NAME}, time=${TIME_LIMIT})"
