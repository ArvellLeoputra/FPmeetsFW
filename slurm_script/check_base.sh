#!/bin/bash

# Check results for base FP-FW variant: manhattan + vanilla + unitary + binary

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"

DEF_GLOBAL_TIME_LIMIT=$(grep 'const DEF_GLOBAL_TIME_LIMIT' "$BASE_DIR/dependencies.jl" | grep -oE '[0-9]+(\.[0-9]+)?')
DEF_SCIP_TIME_LIMIT=$(grep 'const DEF_SCIP_TIME_LIMIT' "$BASE_DIR/dependencies.jl" | grep -oE '[0-9]+(\.[0-9]+)?')

NORM="manhattan"
VARIANT="vanilla"
LS="unitary"

FOLDER="base"
NAME="${NORM}_${VARIANT}_${LS}"

OUTPUT_DIR="$BASE_DIR/$FOLDER/output"
RESULT_DIR="$BASE_DIR/$FOLDER/result"
rm -rf "$RESULT_DIR"
mkdir -p "$RESULT_DIR"

SUMMARY="$RESULT_DIR/solution_summary.txt"
DETAILED="$RESULT_DIR/detailed_results.txt"
FOUND="$RESULT_DIR/solutions_found.txt"
FAILED="$RESULT_DIR/failed_runs.txt"

CSV_HEADER="Instance,BinaryVars,IntegerVars,SolutionFound,TotalTime,FPIterations,FWIterations,Restarts,Objective,Gap,FailureReason,FailureType,ProjectionNorm,FWVariant,LineSearch,RandRound,RandFeasCheck,WarmStart"
echo "$CSV_HEADER" > "$CSV_FILE"

echo "FPFW Solution Analysis [$NAME] - $(date)" > "$SUMMARY"
echo "==========================================================" >> "$SUMMARY"

echo "DETAILED RESULTS" > "$DETAILED"
echo "==========================================================" >> "$DETAILED"
echo "  Projection norm:    $NORM" >> "$DETAILED"
echo "  FW variant:         $VARIANT" >> "$DETAILED"
echo "  Line search:        $LS" >> "$DETAILED"
echo "==========================================================" >> "$DETAILED"
echo "" >> "$DETAILED"

echo "Test cases with solutions found:" > "$FOUND"
echo "=============================================================================================================================" >> "$FOUND"
printf "%-35s %-8s %-8s %-12s %-10s %-10s %-10s %-20s %-10s\n" "Instance" "Bin" "Int" "Time (s)" "FP Iters" "FW Iters" "Restarts" "Objective" "Gap (%)" >> "$FOUND"
echo "=============================================================================================================================" >> "$FOUND"

echo "Failed/Interrupted runs:" > "$FAILED"
echo "===============================================================================================================================" >> "$FAILED"
printf "%-35s %-8s %-8s %-15s %-50s\n" "Instance" "Bin" "Int" "Runtime (s)" "Failure Reason" >> "$FAILED"
echo "===============================================================================================================================" >> "$FAILED"

total_count=0
found_count=0
failed_count=0
scip_timelimit_count=0
global_timelimit_count=0
restart_limit_count=0
fw_infeasible_count=0
other_failure_count=0
rr_found_count=0
found_binary=0
found_ginteger=0
failed_binary=0
failed_ginteger=0

declare -a found_times=()
declare -a found_fw_iters=()
declare -a found_fp_iters=()
declare -a found_restarts=()

for output_file in "$OUTPUT_DIR"/*.out; do
    if [ ! -f "$output_file" ]; then continue; fi

    total_count=$((total_count + 1))
    instance_name=$(basename "$output_file" .out)

    projection_norm=$(grep "projectionNorm =" "$output_file" | head -1 | awk '{print $3}')
    if [ -z "$projection_norm" ]; then projection_norm="$NORM"; fi
    fw_variant=$(grep "fwVariant =" "$output_file" | head -1 | awk '{print $3}')
    line_search=$(grep "lineSearch =" "$output_file" | head -1 | awk '{print $3}')
    rand_round=$(grep "randomizedRounding =" "$output_file" | head -1 | awk '{print $3}')
    rand_feas_check=$(grep "randomizedFeasibilityCheck =" "$output_file" | head -1 | awk '{print $3}')
    warm_start=$(grep "warmStart =" "$output_file" | head -1 | awk '{print $3}')

    binary_vars=$(grep "binaryVars =" "$output_file" | tail -1 | awk '{print $3}')
    integer_vars=$(grep "integerVars =" "$output_file" | tail -1 | awk '{print $3}')
    total_time=$(grep "^totalTime =" "$output_file" | tail -1 | awk '{print $3}' | tr -d 's')
    fw_time=$(grep "fwTime =" "$output_file" | tail -1 | awk '{print $3}' | tr -d 's')
    fp_iterations=$(grep "pumpIterations =" "$output_file" | tail -1 | awk '{print $3}')
    fw_iterations=$(grep "fwIterations =" "$output_file" | tail -1 | awk '{print $3}')
    restarts=$(grep "restartCount =" "$output_file" | tail -1 | awk '{print $3}')
    solution_found=$(grep "solFound =" "$output_file" | tail -1 | awk '{print $3}')
    if [ -z "$solution_found" ]; then solution_found="false"; fi

    exit_reason_str=$(grep "exitReason =" "$output_file" | tail -1 | sed 's/exitReason = //')

    scip_timelimit=""
    global_timelimit=""
    restart_limit=""
    fw_infeasible=""
    rr_found_flag=""

    case "$exit_reason_str" in
        *"SCIP time limit"*)               scip_timelimit="yes" ;;
        *"global time limit"*)             global_timelimit="yes" ;;
        *"cycled"*)                        restart_limit="yes" ;;
        *"outside the feasible polytope"*) fw_infeasible="yes" ;;
        *"randomized rounding"*)           rr_found_flag="yes" ;;
    esac

    if [ "$solution_found" = "true" ]; then
        objective=$(grep "primalBound =" "$output_file" | tail -1 | awk '{print $3}')
        gap_line=$(grep "^gap =" "$output_file" | tail -1)
        if echo "$gap_line" | grep -q "Infinite"; then
            gap="Infinite"
        else
            gap=$(echo "$gap_line" | awk '{print $3}')
        fi
    else
        objective="N/A"
        gap="N/A"
    fi

    failure_reason="None"
    failure_type="None"
    this_was_found=0

    if [ -n "$scip_timelimit" ]; then
        failure_reason="SCIP_TIME_LIMIT (${DEF_SCIP_TIME_LIMIT}s)"
        failure_type="SCIP_TIMELIMIT"
        failed_count=$((failed_count + 1))
        scip_timelimit_count=$((scip_timelimit_count + 1))
        printf "%-35s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${total_time:-N/A}" "SCIP time limit (${DEF_SCIP_TIME_LIMIT}s)" >> "$FAILED"

    elif [ -n "$global_timelimit" ]; then
        failure_reason="GLOBAL_TIME_LIMIT (${DEF_GLOBAL_TIME_LIMIT}s)"
        failure_type="GLOBAL_TIMELIMIT"
        failed_count=$((failed_count + 1))
        global_timelimit_count=$((global_timelimit_count + 1))
        printf "%-35s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${total_time:-N/A}" "Global time limit (${DEF_GLOBAL_TIME_LIMIT}s)" >> "$FAILED"

    elif [ "$solution_found" = "true" ]; then
        found_count=$((found_count + 1))
        this_was_found=1
        [ -n "$rr_found_flag" ] && rr_found_count=$((rr_found_count + 1))
        found_times+=("$total_time")
        found_fw_iters+=("$fw_iterations")
        found_fp_iters+=("$fp_iterations")
        found_restarts+=("$restarts")
        printf "%-35s %-8s %-8s %-12s %-10s %-10s %-10s %-12s %-8s\n" \
            "$instance_name" "$binary_vars" "${integer_vars:-0}" "$total_time" "$fp_iterations" "$fw_iterations" "$restarts" "$objective" "$gap" >> "$FOUND
        "

    elif [ -n "$fw_infeasible" ]; then
        failure_reason="FW_INFEASIBLE"
        failure_type="FW_INFEASIBLE"
        failed_count=$((failed_count + 1))
        fw_infeasible_count=$((fw_infeasible_count + 1))
        printf "%-35s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${total_time:-N/A}" "FW returns infeasible solution" >> "$FAILED"

    elif [ -n "$restart_limit" ]; then
        failure_reason="RESTART_LIMIT"
        failure_type="RESTART_LIMIT"
        failed_count=$((failed_count + 1))
        restart_limit_count=$((restart_limit_count + 1))
        printf "%-35s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${total_time:-N/A}" "Maximum restarts reached" >> "$FAILED"

    else
        failure_reason="UNKNOWN_ERROR"
        failure_type="UNKNOWN"
        failed_count=$((failed_count + 1))
        other_failure_count=$((other_failure_count + 1))
        printf "%-35s %-8s %-8s %-15s %-50s\n" "$instance_name" "$binary_vars" "${integer_vars:-0}" "${total_time:-N/A}" "Unknown error" >> "$FAILED"
    fi

    if [ "${integer_vars:-0}" -gt 0 ] 2>/dev/null; then
        [ "$this_was_found" = "1" ] && found_ginteger=$((found_ginteger + 1)) || failed_ginteger=$((failed_ginteger + 1))
    else
        [ "$this_was_found" = "1" ] && found_binary=$((found_binary + 1)) || failed_binary=$((failed_binary + 1))
    fi

    echo "${instance_name},${binary_vars},${integer_vars},${solution_found},${total_time},${fp_iterations},${fw_iterations},${restarts},${objective},${gap},${failure_reason},${failure_type},${projection_norm},${fw_variant},${line_search},${rand_round},${rand_feas_check},${warm_start}" >> "$CSV_FILE"

    echo "Instance: ${instance_name}" >> "$DETAILED"
    echo "  Binary variables:  ${binary_vars}" >> "$DETAILED"
    echo "  Integer variables: ${integer_vars}" >> "$DETAILED"
    echo "  Solution found:    ${solution_found}" >> "$DETAILED"
    echo "  Total time:        ${total_time}s" >> "$DETAILED"
    echo "  FP iterations:     ${fp_iterations}" >> "$DETAILED"
    echo "  FW iterations:     ${fw_iterations}" >> "$DETAILED"
    echo "  FW time:           ${fw_time}s" >> "$DETAILED"
    echo "  Restarts:          ${restarts}" >> "$DETAILED"
    if [ "$solution_found" = "true" ]; then
        echo "  Objective:         ${objective}" >> "$DETAILED"
        echo "  Gap:               ${gap}%" >> "$DETAILED"
    else
        echo "  Failure reason:    ${failure_reason}" >> "$DETAILED"
    fi
    echo "" >> "$DETAILED"
done

if [ ${#found_times[@]} -gt 0 ]; then
    sum_time=0
    for t in "${found_times[@]}"; do
        sum_time=$(awk -v s="$sum_time" -v v="$t" 'BEGIN {print s + v}')
    done
    avg_time=$(awk -v s="$sum_time" -v n="${#found_times[@]}" 'BEGIN {printf "%.2f", s / n}')

    sum_fw=0
    for fw in "${found_fw_iters[@]}"; do sum_fw=$((sum_fw + fw)); done
    avg_fw=$((sum_fw / ${#found_fw_iters[@]}))

    sum_fp=0
    for fp in "${found_fp_iters[@]}"; do sum_fp=$((sum_fp + fp)); done
    avg_fp=$((sum_fp / ${#found_fp_iters[@]}))

    sum_restarts=0
    for r in "${found_restarts[@]}"; do sum_restarts=$((sum_restarts + r)); done
    avg_restarts=$(awk -v s="$sum_restarts" -v n="${#found_restarts[@]}" 'BEGIN {printf "%.2f", s / n}')
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
echo "  SCIP time limit:    $scip_timelimit_count"
echo "  Global time limit:  $global_timelimit_count"
echo "  Restart limit:      $restart_limit_count"
echo "  FW infeasible:      $fw_infeasible_count"
echo "  Unknown:            $other_failure_count"
if [ ${#found_times[@]} -gt 0 ]; then
    echo ""
    echo "STATISTICS FOR SUCCESSFUL RUNS"
    echo "=========================================================="
    echo "Average time:          ${avg_time}s"
    echo "Average FP iterations: $avg_fp"
    echo "Average FW iterations: $avg_fw"
    echo "Average restarts:      $avg_restarts"
fi
} | tee -a "$SUMMARY"

echo ""
echo "Results saved to $RESULT_DIR"
echo "CSV: $CSV_FILE"
