#!/bin/bash

# Submit base FP-FW variant: euclidean + away + secant
# Sweeps 4 combinations of rand_round and warm_start

PROJECT_DIR="/home/htc/aleoputra/project/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/selection_benchmark"

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR"
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

NORM="euclidean"
VARIANT="away"
LS="secant"

# (rand_round, warm_start) combinations
RR_VALUES=("false" "true"  "false" "true")
WS_VALUES=("false" "false" "true"  "true")

for i in "${!RR_VALUES[@]}"; do
    RR="${RR_VALUES[$i]}"
    WS="${WS_VALUES[$i]}"
    NAME="${NORM}_${VARIANT}_${LS}_rr${RR}_ws${WS}"

    OUT_DIR="$PROJECT_DIR/run_focused/$NAME/output"
    ERR_DIR="$PROJECT_DIR/run_focused/$NAME/error"
    rm -rf "$OUT_DIR" "$ERR_DIR"
    mkdir -p "$OUT_DIR" "$ERR_DIR"

    sbatch <<EOF || { echo "ERROR: sbatch failed for $NAME"; continue; }
#!/bin/bash
#SBATCH --job-name=fpfw_focused_rr${RR}_ws${WS}
#SBATCH --time=20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=big
#SBATCH --constraint=virtual
#SBATCH --array=1-${NUM_INSTANCES}
#SBATCH --output=$OUT_DIR/job_%A_%a.out
#SBATCH --error=$ERR_DIR/job_%A_%a.err

INSTANCE=\$(ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort | sed -n "\${SLURM_ARRAY_TASK_ID}p")

if [ -z "\$INSTANCE" ]; then
    echo "Error: No instance found for task ID \${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INSTANCE_PATH="$INSTANCE_DIR/\$INSTANCE"
echo "Running instance: \$INSTANCE_PATH"
echo "SLURM task ID: \${SLURM_ARRAY_TASK_ID}"

julia --project=$PROJECT_DIR $PROJECT_DIR/run_test.jl "\$INSTANCE_PATH" $NORM $VARIANT $LS $RR $WS
EOF

    echo "Submitted: $NAME"
done

echo "All ${#RR_VALUES[@]} combinations submitted (${NUM_INSTANCES} instances each)."
