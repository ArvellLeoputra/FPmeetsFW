#!/bin/bash

# Submit FP-FW variant: manhattan + vanilla + unitary, fwMaxIterations=1, rfc=true  # UPDATE

PROJECT_DIR="/home/htc/aleoputra/project"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/instances/miplib_selected"  # UPDATE
COMP_RESULT="$PROJECT_DIR/compResult"
CONFIG="$FPFW_DIR/settings/run6.cfg"

NUM_INSTANCES=$(ls "$INSTANCE_DIR" | grep -cE '\.mps(\.gz)?$')
if [ "$NUM_INSTANCES" -eq 0 ]; then
    echo "Error: No .mps/.mps.gz instances found in $INSTANCE_DIR"
    exit 1
fi
echo "Found $NUM_INSTANCES instances in $INSTANCE_DIR"

FOLDER="fpfw_run6"  # UPDATE

OUT_DIR="$COMP_RESULT/$FOLDER/output"
ERR_DIR="$COMP_RESULT/$FOLDER/error"
rm -rf "$OUT_DIR" "$ERR_DIR"
mkdir -p "$OUT_DIR" "$ERR_DIR"

sbatch <<EOF || { echo "ERROR: sbatch failed for $FOLDER"; exit 1; }
#!/bin/bash
#SBATCH --job-name=${FOLDER}
#SBATCH --time=30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --array=1-${NUM_INSTANCES}
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

export PATH="/home/htc/aleoputra/scratch/julia-1.12.6/bin:\$PATH"

INSTANCE=\$(ls "$INSTANCE_DIR" | grep -E '\.mps(\.gz)?$' | sort | sed -n "\${SLURM_ARRAY_TASK_ID}p")

if [ -z "\$INSTANCE" ]; then
    exit 1
fi

BASENAME=\$(basename "\$INSTANCE" .mps.gz)
BASENAME=\$(basename "\$BASENAME" .mps)
exec > "$OUT_DIR/\${BASENAME}.out" 2> "$ERR_DIR/\${BASENAME}.err"

INSTANCE_PATH="$INSTANCE_DIR/\$INSTANCE"
echo "Running instance: \$INSTANCE_PATH"
echo "SLURM task ID: \${SLURM_ARRAY_TASK_ID}"
echo "Config: $CONFIG"

julia --project=$FPFW_DIR $FPFW_DIR/main.jl "\$INSTANCE_PATH" "$CONFIG"
EOF

echo "Submitted: $FOLDER (${NUM_INSTANCES} instances)"
