#!/bin/bash

# Submit all 28 FP-FW experiment combinations (7 norm/linesearch x 4 FW variants)
# Each job runs as an array over 20 instances in selection_benchmark

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"
INSTANCE_DIR="$BASE_DIR/selection_benchmark"
PROJECT_DIR="$BASE_DIR"

NORMS=("manhattan" "euclidean" "abssmooth" "euclidean" "euclidean" "abssmooth" "abssmooth")
LINESEARCHES=("agnostic" "agnostic" "agnostic" "adaptive" "secant" "adaptive" "secant")
VARIANTS=("vanilla" "away" "blended_pairwise" "blended")

# Short name mappings for job names and directories
declare -A NORM_SHORT=( ["manhattan"]="man" ["euclidean"]="euc" ["abssmooth"]="abs" )
declare -A LS_SHORT=( ["agnostic"]="agn" ["adaptive"]="adp" ["secant"]="sec" )
declare -A VAR_SHORT=( ["vanilla"]="van" ["away"]="awy" ["blended_pairwise"]="bpw" ["blended"]="bld" )

for i in "${!NORMS[@]}"; do
    NORM="${NORMS[$i]}"
    LS="${LINESEARCHES[$i]}"

    for VARIANT in "${VARIANTS[@]}"; do
        NAME="${NORM}_${VARIANT}_${LS}"
        SHORT="${NORM_SHORT[$NORM]}_${VAR_SHORT[$VARIANT]}_${LS_SHORT[$LS]}"

        OUT_DIR="$BASE_DIR/run/$NAME/output"
        ERR_DIR="$BASE_DIR/run/$NAME/error"
        mkdir -p "$OUT_DIR" "$ERR_DIR"

        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=fpfw_$SHORT
#SBATCH --time=10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
#SBATCH --array=1-20
#SBATCH --output=$OUT_DIR/job_%A_%a.out
#SBATCH --error=$ERR_DIR/job_%A_%a.err

INSTANCE=\$(ls "$INSTANCE_DIR" | sort | sed -n "\${SLURM_ARRAY_TASK_ID}p")

if [ -z "\$INSTANCE" ]; then
    echo "Error: No instance found for task ID \${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INSTANCE_PATH="$INSTANCE_DIR/\$INSTANCE"
echo "Running instance: \$INSTANCE_PATH"
echo "SLURM task ID: \${SLURM_ARRAY_TASK_ID}"

julia --project $PROJECT_DIR/run_test.jl "\$INSTANCE_PATH" $NORM 0.5 $VARIANT $LS
EOF

        echo "Submitted: $NAME"
    done
done

echo "All 28 jobs submitted."
