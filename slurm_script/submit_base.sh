#!/bin/bash

# Submit base FP-FW variant: manhattan + vanilla + unitary

PROJECT_DIR="/home/htc/aleoputra/project/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/benchmark_full"

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR"
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

NORM="manhattan"
VARIANT="vanilla"
LS="unitary"
NAME="run_base"

OUT_DIR="$PROJECT_DIR/$NAME/output"
ERR_DIR="$PROJECT_DIR/$NAME/error"
rm -rf "$OUT_DIR" "$ERR_DIR"
mkdir -p "$OUT_DIR" "$ERR_DIR"

sbatch <<EOF || { echo "ERROR: sbatch failed for $NAME"; exit 1; }
#!/bin/bash
#SBATCH --job-name=fpfw_${NAME}
#SBATCH --time=20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=big
#SBATCH --constraint=virtual
#SBATCH --array=1-${NUM_INSTANCES}
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

INSTANCE=\$(ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort | sed -n "\${SLURM_ARRAY_TASK_ID}p")

if [ -z "\$INSTANCE" ]; then
    exit 1
fi

BASENAME=\$(basename "\$INSTANCE" .mps)
exec > "$OUT_DIR/\${BASENAME}.out" 2> "$ERR_DIR/\${BASENAME}.err"

INSTANCE_PATH="$INSTANCE_DIR/\$INSTANCE"
echo "Running instance: \$INSTANCE_PATH"
echo "SLURM task ID: \${SLURM_ARRAY_TASK_ID}"

julia --project=$PROJECT_DIR $PROJECT_DIR/run_test.jl "\$INSTANCE_PATH" $NORM $VARIANT $LS
EOF

echo "Submitted: $NAME (${NUM_INSTANCES} instances)"