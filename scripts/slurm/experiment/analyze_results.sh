#!/bin/bash

# Check results for a FP-FW variant. Usage: ./analyze_results.sh [cfgName] [folder]
# cfgName defaults to "fpfw" (settings/fpfw.cfg), folder defaults to "fpfw_baseline"

PROJECT_DIR="${PROJECT_DIR:-/home/htc/aleoputra/project}"
FPFW_DIR="$PROJECT_DIR/FPmeetsFW"
COMP_RESULT="$PROJECT_DIR/compResult"
SETTINGS_DIR="$FPFW_DIR/settings"

CFG_NAME="${1:-fpfw}"
FOLDER="${2:-fpfw_baseline}"

RESULT_DIR="$COMP_RESULT/$FOLDER"

if ! command -v jq > /dev/null 2>&1; then
    echo "Error: jq not found; required to parse results.json" >&2
    exit 1
fi

# Prefer the config archived at submission time so results and settings stay
# in sync even if settings/${CFG_NAME}.cfg has since been edited. Fall back
# to the live settings file for runs submitted before this archiving existed.
RUN_CONFIG="$RESULT_DIR/config.cfg"
if [ -f "$RUN_CONFIG" ]; then
    CONFIG="$RUN_CONFIG"
else
    CONFIG="$SETTINGS_DIR/${CFG_NAME}.cfg"
fi

RUN_NAME=$(grep '^runName=' "$CONFIG" | cut -d'=' -f2)
SEED=$(grep '^seed=' "$CONFIG" | cut -d'=' -f2)
NORM=$(grep '^norm=' "$CONFIG" | cut -d'=' -f2)
VARIANT=$(grep '^fwVariant=' "$CONFIG" | cut -d'=' -f2)
FW_MAX_ITER=$(grep '^fwMaxIterations=' "$CONFIG" | cut -d'=' -f2)
STEP_SIZE=$(grep '^fwStepSize=' "$CONFIG" | cut -d'=' -f2)
TIME_LIMIT=$(grep '^timeLimit=' "$CONFIG" | cut -d'=' -f2)
RAND_ROUND=$(grep '^randomizedRounding=' "$CONFIG" | cut -d'=' -f2)
RAND_FEAS_CHECK=$(grep '^randomFeasibilityCheck=' "$CONFIG" | cut -d'=' -f2)
FW_WARM_START=$(grep '^fwWarmStart=' "$CONFIG" | cut -d'=' -f2)
LMO_WARM_START=$(grep '^lmoWarmStart=' "$CONFIG" | cut -d'=' -f2)
USE_DIVE=$(grep '^useDive=' "$CONFIG" | cut -d'=' -f2)
PRESOLVE=$(grep '^presolve=' "$CONFIG" | cut -d'=' -f2)
VERBOSE=$(grep '^verbose=' "$CONFIG" | cut -d'=' -f2)

RES_DIR="$COMP_RESULT/config_comparison/$FOLDER"
rm -rf "$RES_DIR"
mkdir -p "$RES_DIR"

SUMMARY="$RES_DIR/solution_summary.txt"
DETAILED="$RES_DIR/detailed_results.txt"
FOUND="$RES_DIR/solutions_found.txt"
FAILED="$RES_DIR/failed_runs.txt"

echo "FPFW Solution Analysis [$FOLDER] - $(date)" > "$SUMMARY"
echo "==========================================================" >> "$SUMMARY"

echo "DETAILED RESULTS" > "$DETAILED"
echo "==========================================================" >> "$DETAILED"
echo "  Run name:           $RUN_NAME" >> "$DETAILED"
echo "  Seed:               $SEED" >> "$DETAILED"
echo "  Projection norm:    $NORM" >> "$DETAILED"
echo "  FW variant:         $VARIANT" >> "$DETAILED"
echo "  FW max iterations:  $FW_MAX_ITER" >> "$DETAILED"
echo "  Step size:          $STEP_SIZE" >> "$DETAILED"
echo "  Time limit:         ${TIME_LIMIT}s" >> "$DETAILED"
echo "  Rand rounding:      $RAND_ROUND" >> "$DETAILED"
echo "  Rand feas check:    $RAND_FEAS_CHECK" >> "$DETAILED"
echo "  FW warm start:      $FW_WARM_START" >> "$DETAILED"
echo "  LMO warm start:     $LMO_WARM_START" >> "$DETAILED"
echo "  Use dive:           $USE_DIVE" >> "$DETAILED"
echo "  Presolve:           $PRESOLVE" >> "$DETAILED"
echo "  Verbose:            $VERBOSE" >> "$DETAILED"
echo "==========================================================" >> "$DETAILED"
echo "" >> "$DETAILED"

echo "Test cases with solutions found:" > "$FOUND"
echo "=====================================================================================================================================" >> "$FOUND"
printf "%-35s %-8s %-8s %-8s %-12s %-12s %-10s %-10s %-10s %-10s %-16s %-10s\n" "Instance" "Bin" "Int" "Cont" "Total(s)" "HeurT(s)" "FP Iters" "FW Iters" "Restarts" "Perturbs" "Objective" "Gap (%)" >> "$FOUND"
echo "=====================================================================================================================================" >> "$FOUND"

echo "Failed/Interrupted runs:" > "$FAILED"
echo "===============================================================================================================================" >> "$FAILED"
printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "Instance" "Bin" "Int" "Cont" "heurTime (s)" "Failure Reason" >> "$FAILED"
echo "===============================================================================================================================" >> "$FAILED"

fmt_obj() {
    awk -v v="$1" 'BEGIN {
        x = v + 0
        if (x == int(x)) printf "%d", x
        else printf "%.4g", x
    }'
}

total_count=0
found_count=0
failed_count=0
scip_timelimit_count=0
timelimit_count=0
iterlimit_count=0
fw_infeasible_count=0
other_failure_count=0
rr_found_count=0
found_binary=0
found_ginteger=0
failed_binary=0
failed_ginteger=0

declare -a found_times=()
declare -a found_heur_times=()
declare -a found_fw_iters=()
declare -a found_fp_iters=()
declare -a found_restarts=()
declare -a found_perturbations=()

for instance_dir in "$RESULT_DIR"/*/; do
    instance_dir="${instance_dir%/}"
    instance_name=$(basename "$instance_dir")
    [ "$instance_name" = "result" ] && continue

    output_file="$instance_dir/slurm_job.out"
    results_json="$instance_dir/results.json"
    if [ ! -f "$output_file" ]; then continue; fi

    total_count=$((total_count + 1))

    binary_vars=$(grep "binaryVars =" "$output_file" | tail -1 | awk '{print $3}')
    integer_vars=$(grep "integerVars =" "$output_file" | tail -1 | awk '{print $3}')
    continuous_vars=$(grep "continuousVars =" "$output_file" | tail -1 | awk '{print $3}')

    if [ -f "$results_json" ]; then
        exit_reason_code=$(jq -r '.exitReason' "$results_json")
        solution_found=$(jq -r '.solutionFound' "$results_json")
        total_time=$(awk -v t="$(jq -r '.totalTime' "$results_json")" 'BEGIN { printf "%.2f", t }')
        heur_time=$(awk -v t="$(jq -r '.heurTime' "$results_json")" 'BEGIN { printf "%.2f", t }')
        fw_time=$(awk -v t="$(jq -r '.fwTime' "$results_json")" 'BEGIN { printf "%.2f", t }')
        fp_iterations=$(jq -r '.pumpIterations' "$results_json")
        fw_iterations=$(jq -r '.fwIterations' "$results_json")
        restarts=$(jq -r '.restartCount' "$results_json")
        perturbations=$(jq -r '.perturbCount' "$results_json")
        objective=$(jq -r '.primalBound' "$results_json")
        gap_frac=$(jq -r '.gap' "$results_json")
        gap=$(awk -v g="$gap_frac" 'BEGIN { printf "%.2f", g * 100 }')
        # heurTime/fwTime/etc. are meaningless when SCIP finished before the
        # heuristic ever ran; blank them out so they render as N/A below.
        case "$exit_reason_code" in
            SCIP_*) heur_time="" ;;
        esac
    else
        # No results.json means the run crashed before writing it.
        exit_reason_code="UNKNOWN"
        solution_found="false"
        total_time=""
        heur_time=""
        fw_time=""
        fp_iterations=""
        fw_iterations=""
        restarts=""
        perturbations=""
        objective="N/A"
        gap="N/A"
    fi

    if [ "$solution_found" != "true" ]; then
        objective="N/A"
        gap="N/A"
    fi

    failure_reason="None"
    failure_type="None"
    this_was_found=0

    if [ "$exit_reason_code" = "SCIP_TIME_LIMIT" ]; then
        failure_reason="SCIP_TIME_LIMIT (${TIME_LIMIT}s)"
        failure_type="SCIP_TIMELIMIT"
        failed_count=$((failed_count + 1))
        scip_timelimit_count=$((scip_timelimit_count + 1))
        printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "${heur_time:-N/A}" "SCIP time limit (${TIME_LIMIT}s)" >> "$FAILED"

    elif [ "$exit_reason_code" = "TIME_LIMIT" ]; then
        failure_reason="TIME_LIMIT (${TIME_LIMIT}s)"
        failure_type="TIMELIMIT"
        failed_count=$((failed_count + 1))
        timelimit_count=$((timelimit_count + 1))
        printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "${heur_time:-N/A}" "Time limit (${TIME_LIMIT}s)" >> "$FAILED"

    elif [ "$exit_reason_code" = "ITER_LIMIT" ]; then
        failure_reason="ITER_LIMIT"
        failure_type="ITERLIMIT"
        failed_count=$((failed_count + 1))
        iterlimit_count=$((iterlimit_count + 1))
        printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "${heur_time:-N/A}" "Iteration limit reached" >> "$FAILED"

    elif [ "$solution_found" = "true" ]; then
        found_count=$((found_count + 1))
        this_was_found=1
        [ "$exit_reason_code" = "SOLUTION_RR" ] && rr_found_count=$((rr_found_count + 1))
        found_times+=("$total_time")
        found_heur_times+=("${heur_time:-0}")
        found_fw_iters+=("$fw_iterations")
        found_fp_iters+=("$fp_iterations")
        found_restarts+=("$restarts")
        found_perturbations+=("$perturbations")
        fmt_objective=$(fmt_obj "$objective")
        printf "%-35s %-8s %-8s %-8s %-12s %-12s %-10s %-10s %-10s %-10s %-16s %-10s\n" \
            "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "$total_time" "${heur_time:-N/A}" "$fp_iterations" "$fw_iterations" "$restarts" "$perturbations" "$fmt_objective" "$gap" >> "$FOUND"

    elif [ "$exit_reason_code" = "INFEASIBLE_FW" ]; then
        failure_reason="FW_INFEASIBLE"
        failure_type="FW_INFEASIBLE"
        failed_count=$((failed_count + 1))
        fw_infeasible_count=$((fw_infeasible_count + 1))
        printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "${heur_time:-N/A}" "FW returns infeasible solution" >> "$FAILED"

    else
        failure_reason="UNKNOWN_ERROR"
        failure_type="UNKNOWN"
        failed_count=$((failed_count + 1))
        other_failure_count=$((other_failure_count + 1))
        printf "%-35s %-8s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${continuous_vars:-0}" "${heur_time:-N/A}" "${exit_reason_code}" >> "$FAILED"
    fi

    if [ "${integer_vars:-0}" -gt 0 ] 2>/dev/null; then
        [ "$this_was_found" = "1" ] && found_ginteger=$((found_ginteger + 1)) || failed_ginteger=$((failed_ginteger + 1))
    else
        [ "$this_was_found" = "1" ] && found_binary=$((found_binary + 1)) || failed_binary=$((failed_binary + 1))
    fi

    echo "Instance: ${instance_name}" >> "$DETAILED"
    printf "  %-18s%s\n" "Binary Vars:" "${binary_vars}" >> "$DETAILED"
    printf "  %-18s%s\n" "Integer Vars:" "${integer_vars}" >> "$DETAILED"
    printf "  %-18s%s\n" "Continuous Vars:" "${continuous_vars}" >> "$DETAILED"
    printf "  %-18s%s\n" "Solution found:" "${solution_found}" >> "$DETAILED"
    printf "  %-18s%s\n" "Total time:" "${total_time}s" >> "$DETAILED"
    printf "  %-18s%s\n" "Heur time:" "${heur_time:-N/A}s" >> "$DETAILED"
    printf "  %-18s%s\n" "FP iterations:" "${fp_iterations}" >> "$DETAILED"
    printf "  %-18s%s\n" "FW iterations:" "${fw_iterations}" >> "$DETAILED"
    printf "  %-18s%s\n" "FW time:" "${fw_time}s" >> "$DETAILED"
    printf "  %-18s%s\n" "Perturbations:" "${perturbations}" >> "$DETAILED"
    printf "  %-18s%s\n" "Restarts:" "${restarts}" >> "$DETAILED"
    printf "  %-18s%s\n" "Exit reason:" "${exit_reason_code}" >> "$DETAILED"
    if [ "$solution_found" = "true" ]; then
        printf "  %-18s%s\n" "Objective:" "${objective}" >> "$DETAILED"
        printf "  %-18s%s%%\n" "Gap:" "${gap}" >> "$DETAILED"
    else
        printf "  %-18s%s\n" "Gap:" "N/A" >> "$DETAILED"
    fi
    echo "" >> "$DETAILED"
done

if [ ${#found_times[@]} -gt 0 ]; then
    sum_time=0
    for t in "${found_times[@]}"; do
        sum_time=$(awk -v s="$sum_time" -v v="$t" 'BEGIN {print s + v}')
    done
    avg_time=$(awk -v s="$sum_time" -v n="${#found_times[@]}" 'BEGIN {printf "%.2f", s / n}')

    sum_heur_time=0
    for t in "${found_heur_times[@]}"; do
        sum_heur_time=$(awk -v s="$sum_heur_time" -v v="$t" 'BEGIN {print s + v}')
    done
    avg_heur_time=$(awk -v s="$sum_heur_time" -v n="${#found_heur_times[@]}" 'BEGIN {printf "%.2f", s / n}')

    sum_fw=0
    for fw in "${found_fw_iters[@]}"; do sum_fw=$((sum_fw + fw)); done
    avg_fw=$((sum_fw / ${#found_fw_iters[@]}))

    sum_fp=0
    for fp in "${found_fp_iters[@]}"; do sum_fp=$((sum_fp + fp)); done
    avg_fp=$((sum_fp / ${#found_fp_iters[@]}))

    sum_restarts=0
    for r in "${found_restarts[@]}"; do sum_restarts=$((sum_restarts + r)); done
    avg_restarts=$(awk -v s="$sum_restarts" -v n="${#found_restarts[@]}" 'BEGIN {printf "%.2f", s / n}')

    sum_perturbations=0
    for p in "${found_perturbations[@]}"; do sum_perturbations=$((sum_perturbations + p)); done
    avg_perturbations=$(awk -v s="$sum_perturbations" -v n="${#found_perturbations[@]}" 'BEGIN {printf "%.2f", s / n}')
fi

{
echo ""
echo "OVERALL STATISTICS"
echo "=========================================================="
echo "Total instances:    $total_count"
echo "Solutions found:    $found_count ($(awk -v f="$found_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)"
echo "  Via rand round:   $rr_found_count"
echo "  Via FPFW loop:    $((found_count - rr_found_count))"
echo "Failed:             $failed_count ($(awk -v f="$failed_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)"
echo ""
echo "BREAKDOWN BY INSTANCE TYPE"
echo "=========================================================="
echo "Binary-only:     found=$found_binary  failed=$failed_binary"
echo "General integer: found=$found_ginteger  failed=$failed_ginteger"
echo ""
echo "FAILURE BREAKDOWN"
echo "=========================================================="
echo "  Time limit:         $timelimit_count"
echo "  SCIP time limit:    $scip_timelimit_count"
echo "  Iteration limit:    $iterlimit_count"
echo "  FW infeasible:      $fw_infeasible_count"
echo "  Unknown:            $other_failure_count"
if [ ${#found_times[@]} -gt 0 ]; then
    echo ""
    echo "STATISTICS FOR SUCCESSFUL RUNS"
    echo "=========================================================="
    echo "Average total time:    ${avg_time}s"
    echo "Average heur time:     ${avg_heur_time}s"
    echo "Average FP iterations: $avg_fp"
    echo "Average FW iterations: $avg_fw"
    echo "Average perturbations: $avg_perturbations"
    echo "Average restarts:      $avg_restarts"
fi
} | tee -a "$SUMMARY"

echo ""
echo "Results saved to $RES_DIR"
