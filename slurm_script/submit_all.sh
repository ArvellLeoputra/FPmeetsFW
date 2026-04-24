#!/bin/bash

# Submit all FP-FW experiment combinations
# Each job runs as an array over all instances in selection_benchmark

PROJECT_DIR="/home/htc/aleoputra/project/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/selection_benchmark"

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR"
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

NORMS=("manhattan" "euclidean" "abssmooth" "euclidean" "euclidean" "abssmooth" "abssmooth")
LINESEARCHES=("agnostic" "agnostic" "agnostic" "adaptive" "secant" "adaptive" "secant")
VARIANTS=("vanilla" "away" "blended_pairwise" "blended")
PRESOLVES=("false" "true")

for PRESOLVE in "${PRESOLVES[@]}"; do
    for i in "${!NORMS[@]}"; do
        NORM="${NORMS[$i]}"
        LS="${LINESEARCHES[$i]}"

        for VARIANT in "${VARIANTS[@]}"; do
            # blended (BCG) requires curvature-aware line search; skip incompatible pairs
            if [ "$VARIANT" = "blended" ] && { [ "$LS" = "agnostic" ] || [ "$LS" = "secant" ]; }; then
                continue
            fi

            NAME="${NORM}_${VARIANT}_${LS}_presolve_${PRESOLVE}"

            OUT_DIR="$PROJECT_DIR/run/$NAME/output"
            ERR_DIR="$PROJECT_DIR/run/$NAME/error"
            rm -rf "$OUT_DIR" "$ERR_DIR"
            mkdir -p "$OUT_DIR" "$ERR_DIR"

            sbatch <<EOF || { echo "ERROR: sbatch failed for $NAME"; continue; }
#!/bin/bash
#SBATCH --job-name=fpfw_${NORM}_${VARIANT}_${LS}_${PRESOLVE}
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

julia --project=$PROJECT_DIR $PROJECT_DIR/run_test.jl "\$INSTANCE_PATH" $NORM 0.5 $VARIANT $LS $PRESOLVE
EOF

            echo "Submitted: $NAME"
        done
    done
done

echo "All $((${#NORMS[@]} * ${#VARIANTS[@]} * ${#PRESOLVES[@]})) jobs submitted (${NUM_INSTANCES} instances each)."
