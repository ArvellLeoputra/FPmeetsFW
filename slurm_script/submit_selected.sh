#!/bin/bash

# Submit selected FP-FW combinations
# Total: 8 combinations

PROJECT_DIR="/home/htc/aleoputra/project/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/selection_benchmark"

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR"
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

# Explicit list of (NORM, VARIANT, LS) combinations
NORMS=("manhattan" "manhattan" "euclidean" "euclidean" "euclidean" "smooth_manhattan" "smooth_manhattan" "smooth_manhattan")
VARIANTS=("away" "blended_pairwise" "away" "away" "blended" "away" "away" "blended")
LINESEARCHES=("agnostic" "agnostic" "secant" "adaptive" "adaptive" "secant" "adaptive" "adaptive")

for i in "${!NORMS[@]}"; do
    NORM="${NORMS[$i]}"
    VARIANT="${VARIANTS[$i]}"
    LS="${LINESEARCHES[$i]}"
    NAME="${NORM}_${VARIANT}_${LS}"

    OUT_DIR="$PROJECT_DIR/run_selected/$NAME/output"
    ERR_DIR="$PROJECT_DIR/run_selected/$NAME/error"
    rm -rf "$OUT_DIR" "$ERR_DIR"
    mkdir -p "$OUT_DIR" "$ERR_DIR"

    sbatch <<EOF || { echo "ERROR: sbatch failed for $NAME"; continue; }
#!/bin/bash
#SBATCH --job-name=fpfw_${NORM}_${VARIANT}_${LS}
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

julia --project=$PROJECT_DIR $PROJECT_DIR/run_test.jl "\$INSTANCE_PATH" $NORM $VARIANT $LS
EOF

    echo "Submitted: $NAME"
done

echo "All ${#NORMS[@]} combinations submitted (${NUM_INSTANCES} instances each)."
