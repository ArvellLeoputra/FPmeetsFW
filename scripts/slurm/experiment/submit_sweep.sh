#!/bin/bash
#
# submit_sweep.sh - launch the whole {config x seed x instance} comparison as ONE
# interleaved SLURM array job. Because every cell runs in the same job, in the
# same time window, on the same pool of nodes, cross-config and cross-seed
# comparisons are not distorted by "which day / how busy the cluster was".
# (That artifact is exactly why the Aug-24 vs Aug-25 runs weren't comparable.)
#
# Each task writes to:
#   compResult/<folder>_s<seed>/<instance>/{slurm_job.out,slurm_job.err,results.json}
# Bookkeeping for the sweep goes to:
#   compResult/<name>/{sweep_list.tsv,job_script.sh,MANIFEST,JOBID}
#
# Aggregate afterwards with analyze_configs.sh, which globs "<folder>_s*".
#
# Usage:
#   ./submit_sweep.sh [options] <cfg> [<cfg> ...]
#
#   <cfg>   settings/<cfg>.cfg basename. Result folder defaults to "fpfw_<cfg>",
#           except "fpfw" -> "fpfw_baseline". Override per config as "cfg:folder".
#
# Options:
#   -s "LIST"   space-separated seed list        (default: "1 2 3 4 5")
#   -t TIME     SBATCH --time per task           (default: 1:00:00)
#   -j N        max concurrent array tasks (%N)  (default: 64)
#   -n NAME     sweep name / bookkeeping dir     (default: sweep)
#   -S PATH     run each task as "julia --sysimage PATH ..." to skip JIT/load
#               (build it with: julia scripts/sysimage/build.jl)
#   -x          add "#SBATCH --exclusive" (one task per node; use with small -j
#               and a short instance list for a trustworthy-absolute-time run)
#   -y          don't prompt before wiping existing <folder>_s<seed> dirs
#   -d          dry run: build sweep_list.tsv + job_script.sh, don't sbatch
#
# Example - the 6 configs, 5 seeds, 128 wide, with a prebuilt sysimage:
#   julia scripts/sysimage/build.jl
#   ./submit_sweep.sh -j 128 -S ../../../fpfw_sysimage.so fpfw run1 run2 run3 run4 run5
#   # then, when the job finishes:
#   ./analyze_configs.sh fpfw_baseline fpfw_run1 fpfw_run2 fpfw_run3 fpfw_run4 fpfw_run5

set -euo pipefail

usage() {
    sed -n '2,39p' "$0" | sed 's/^#\{1,\} \{0,1\}//; s/^#$//'
}

SEEDS="1 2 3 4 5"
TIME_LIMIT="1:00:00"
THROTTLE=64
SWEEP_NAME="sweep"
EXCLUSIVE=0
ASSUME_YES=0
DRY_RUN=0
SYSIMAGE=""

while getopts "s:t:j:n:S:xyd" opt; do
    case "$opt" in
        s) SEEDS="$OPTARG" ;;
        t) TIME_LIMIT="$OPTARG" ;;
        j) THROTTLE="$OPTARG" ;;
        n) SWEEP_NAME="$OPTARG" ;;
        S) SYSIMAGE="$OPTARG" ;;
        x) EXCLUSIVE=1 ;;
        y) ASSUME_YES=1 ;;
        d) DRY_RUN=1 ;;
        *) usage >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ -n "$SYSIMAGE" ]; then
    case "$SYSIMAGE" in /*) : ;; *) SYSIMAGE="$PWD/$SYSIMAGE" ;; esac
    [ -f "$SYSIMAGE" ] || { echo "Error: sysimage not found: $SYSIMAGE" >&2; exit 1; }
fi

if [ $# -lt 1 ]; then
    echo "Error: need at least one <cfg>" >&2
    usage >&2
    exit 1
fi

case "$THROTTLE" in ''|*[!0-9]*) echo "Error: -j must be a positive integer, got '$THROTTLE'" >&2; exit 1 ;; esac
[ "$THROTTLE" -ge 1 ] || { echo "Error: -j must be >= 1" >&2; exit 1; }

PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
INSTANCE_DIR="$PROJECT_DIR/instances/miplib_selected"
COMP_RESULT="$PROJECT_DIR/compResult"
SETTINGS_DIR="$FPFW_DIR/settings"
JULIA_BIN="${JULIA_BIN:-/home/htc/aleoputra/scratch/julia-1.12.6/bin}"

if [ "$DRY_RUN" != 1 ] && ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: sbatch not found (use -d for a dry run)" >&2
    exit 1
fi

# ---- resolve each <cfg> to (cfgName, folder) and validate the settings file ----
declare -a CFGS=() FOLDERS=()
for a in "$@"; do
    if [[ "$a" == *:* ]]; then
        cfg="${a%%:*}"; folder="${a##*:}"
    elif [ "$a" = "fpfw" ]; then
        cfg="fpfw"; folder="fpfw_baseline"
    else
        cfg="$a"; folder="fpfw_$a"
    fi
    if [ ! -f "$SETTINGS_DIR/$cfg.cfg" ]; then
        echo "Error: config not found: $SETTINGS_DIR/$cfg.cfg" >&2
        exit 1
    fi
    CFGS+=("$cfg"); FOLDERS+=("$folder")
done

# ---- seeds ----
read -ra SEED_ARR <<< "$SEEDS"
[ ${#SEED_ARR[@]} -ge 1 ] || { echo "Error: empty seed list (-s)" >&2; exit 1; }
for s in "${SEED_ARR[@]}"; do
    case "$s" in ''|*[!0-9-]*) echo "Error: seed not an integer: '$s'" >&2; exit 1 ;; esac
done

# ---- instances (sorted, stable order -> stable task ids) ----
mapfile -t INSTANCES < <(ls "$INSTANCE_DIR" 2>/dev/null | grep -E '\.mps(\.gz)?$' | sort || true)
NUM_INST=${#INSTANCES[@]}
[ "$NUM_INST" -gt 0 ] || { echo "Error: no .mps/.mps.gz instances in $INSTANCE_DIR" >&2; exit 1; }

NUM_TASKS=$(( ${#CFGS[@]} * ${#SEED_ARR[@]} * NUM_INST ))

# ---- summary + wipe confirmation ----
declare -a TARGETS=()
for folder in "${FOLDERS[@]}"; do
    for seed in "${SEED_ARR[@]}"; do
        TARGETS+=("${folder}_s${seed}")
    done
done

echo "Sweep '$SWEEP_NAME'"
echo "  configs   : ${CFGS[*]}"
echo "  folders   : ${FOLDERS[*]}"
echo "  seeds     : ${SEED_ARR[*]}"
echo "  instances : $NUM_INST  ($INSTANCE_DIR)"
echo "  tasks     : $NUM_TASKS   (array 1-$NUM_TASKS%$THROTTLE, --time=$TIME_LIMIT/task$([ "$EXCLUSIVE" = 1 ] && echo ', --exclusive'))"
echo "  sysimage  : ${SYSIMAGE:-<none, plain julia>}"
echo "  writes    : $COMP_RESULT/{$(IFS=,; echo "${TARGETS[*]}")}/"

existing=0
for t in "${TARGETS[@]}"; do [ -d "$COMP_RESULT/$t" ] && existing=$((existing + 1)); done
if [ "$existing" -gt 0 ] && [ "$ASSUME_YES" != 1 ] && [ "$DRY_RUN" != 1 ]; then
    if [ -t 0 ]; then
        read -r -p "  $existing of these result dirs exist and will be DELETED. Proceed? [y/N] " ans || ans=""
        case "$ans" in [yY]|[yY][eE][sS]) ;; *) echo "Aborted."; exit 1 ;; esac
    else
        echo "Error: $existing target dirs already exist; pass -y to overwrite." >&2
        exit 1
    fi
fi

# ---- build per-(folder,seed) result trees + archived configs, and the task list ----
SWEEP_DIR="$COMP_RESULT/$SWEEP_NAME"
mkdir -p "$SWEEP_DIR"
SWEEP_LIST="$SWEEP_DIR/sweep_list.tsv"
: > "$SWEEP_LIST"

tid=0
for i in "${!CFGS[@]}"; do
    cfg="${CFGS[$i]}"
    folder="${FOLDERS[$i]}"
    src="$SETTINGS_DIR/$cfg.cfg"
    for seed in "${SEED_ARR[@]}"; do
        rd="$COMP_RESULT/${folder}_s${seed}"
        if [ "$DRY_RUN" != 1 ]; then
            rm -rf "$rd"
            mkdir -p "$rd"
            cp "$src" "$rd/config.cfg"
            if grep -q '^seed=' "$rd/config.cfg"; then
                sed -i "s/^seed=.*/seed=$seed/" "$rd/config.cfg"
            else
                echo "seed=$seed" >> "$rd/config.cfg"
            fi
        fi
        for inst in "${INSTANCES[@]}"; do
            tid=$((tid + 1))
            printf '%d\t%s\t%s\t%s\t%s\t%s\n' \
                "$tid" "$folder" "$seed" "$rd" "$inst" "$INSTANCE_DIR/$inst" >> "$SWEEP_LIST"
        done
    done
done

if [ "$tid" -ne "$NUM_TASKS" ]; then
    echo "Error: built $tid task rows but expected $NUM_TASKS" >&2
    exit 1
fi

# ---- render the array job script (literal template + token substitution) ----
JOB="$SWEEP_DIR/job_script.sh"
cat > "$JOB" <<'TEMPLATE'
#!/bin/bash
#SBATCH --job-name=@@NAME@@
#SBATCH --time=@@WALLTIME@@
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=big
#SBATCH --constraint=Gold6338
@@EXCLUSIVE@@
#SBATCH --array=1-@@NTASKS@@%@@THROTTLE@@
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

export PATH="@@JULIA_BIN@@:$PATH"

SWEEP_LIST="@@SWEEP_LIST@@"
FPFW_DIR="@@FPFW_DIR@@"

# task id -> "id  folder  seed  resultDir  instance  instancePath"
LINE=$(awk -F'\t' -v id="${SLURM_ARRAY_TASK_ID}" '$1 == id {print; exit}' "$SWEEP_LIST")
if [ -z "$LINE" ]; then
    echo "No sweep_list row for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
IFS=$'\t' read -r _ FOLDER SEED RESULT_DIR INSTANCE INSTANCE_PATH <<< "$LINE"

BASENAME=$(basename "$INSTANCE" .mps.gz)
BASENAME=$(basename "$BASENAME" .mps)
INSTANCE_RESULT_DIR="$RESULT_DIR/$BASENAME"
mkdir -p "$INSTANCE_RESULT_DIR"

exec > "$INSTANCE_RESULT_DIR/slurm_job.out" 2> "$INSTANCE_RESULT_DIR/slurm_job.err"

echo "Sweep task:  ${SLURM_ARRAY_TASK_ID}"
echo "Config:      $FOLDER (seed $SEED)"
echo "Config file: $RESULT_DIR/config.cfg"
echo "Instance:    $INSTANCE_PATH"
echo "Timestamp:   $(date '+%Y-%m-%d %H:%M:%S')"

julia @@SYSIMAGE@@--project="$FPFW_DIR" "$FPFW_DIR/main.jl" \
    "$INSTANCE_PATH" "$RESULT_DIR/config.cfg" "resultsDir=$INSTANCE_RESULT_DIR"
TEMPLATE

if [ "$EXCLUSIVE" = 1 ]; then
    excl_line='#SBATCH --exclusive'
else
    excl_line='# (shared nodes; pass -x to submit_sweep.sh for --exclusive)'
fi

# "--sysimage <path> " (trailing space) or "" -> "julia [--sysimage ...]--project=..."
if [ -n "$SYSIMAGE" ]; then
    sysimage_tok="--sysimage $SYSIMAGE "
else
    sysimage_tok=""
fi

# paths here contain no '|', safe as sed delimiter
sed -i \
    -e "s|@@NAME@@|$SWEEP_NAME|g" \
    -e "s|@@WALLTIME@@|$TIME_LIMIT|g" \
    -e "s|@@NTASKS@@|$NUM_TASKS|g" \
    -e "s|@@THROTTLE@@|$THROTTLE|g" \
    -e "s|@@EXCLUSIVE@@|$excl_line|g" \
    -e "s|@@JULIA_BIN@@|$JULIA_BIN|g" \
    -e "s|@@SWEEP_LIST@@|$SWEEP_LIST|g" \
    -e "s|@@FPFW_DIR@@|$FPFW_DIR|g" \
    -e "s|@@SYSIMAGE@@|$sysimage_tok|g" \
    "$JOB"
chmod +x "$JOB"

# ---- manifest ----
{
    echo "sweep_name : $SWEEP_NAME"
    echo "created    : $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "configs    : ${CFGS[*]}"
    echo "folders    : ${FOLDERS[*]}"
    echo "seeds      : ${SEED_ARR[*]}"
    echo "instances  : $NUM_INST ($INSTANCE_DIR)"
    echo "tasks      : $NUM_TASKS"
    echo "throttle   : $THROTTLE"
    echo "walltime   : $TIME_LIMIT"
    echo "exclusive  : $EXCLUSIVE"
    echo "sysimage   : ${SYSIMAGE:-<none>}"
    echo "sweep_list : $SWEEP_LIST"
    echo "job_script : $JOB"
} > "$SWEEP_DIR/MANIFEST"

ANALYZE_CMD="./analyze_configs.sh ${FOLDERS[*]}"

if [ "$DRY_RUN" = 1 ]; then
    echo
    echo "[dry run] built (no sbatch, result dirs untouched):"
    echo "  $SWEEP_LIST   ($NUM_TASKS rows)"
    echo "  $JOB"
    echo "  $SWEEP_DIR/MANIFEST"
    echo "First / last task rows:"
    head -n1 "$SWEEP_LIST" | sed 's/^/  /'
    tail -n1 "$SWEEP_LIST" | sed 's/^/  /'
    echo "Re-run without -d to submit."
    exit 0
fi

JOBID=$(sbatch --parsable "$JOB")
echo "$JOBID" > "$SWEEP_DIR/JOBID"
echo "jobid      : $JOBID" >> "$SWEEP_DIR/MANIFEST"

echo
echo "Submitted array job $JOBID  ($NUM_TASKS tasks, up to $THROTTLE concurrent)"
echo "  watch:      squeue -j $JOBID"
echo "  aggregate:  $ANALYZE_CMD"
