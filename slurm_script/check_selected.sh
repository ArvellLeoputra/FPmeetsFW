#!/bin/bash

# Check results for selected FP-FW combinations
# Generates per-combination result files, combined comparison.csv, and benchmark_summary.txt

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"

DEF_GLOBAL_TIME_LIMIT=$(grep 'const DEF_GLOBAL_TIME_LIMIT' "$BASE_DIR/dependencies.jl" | grep -oE '[0-9]+(\.[0-9]+)?')
DEF_SCIP_TIME_LIMIT=$(grep 'const DEF_SCIP_TIME_LIMIT' "$BASE_DIR/dependencies.jl" | grep -oE '[0-9]+(\.[0-9]+)?')

PRESOLVE="true"

# Explicit list of (NORM, VARIANT, LS) combinations — must match submit_selected.sh
NORMS=("manhattan" "manhattan" "euclidean" "euclidean" "euclidean" "abssmooth" "abssmooth" "abssmooth")
VARIANTS=("away" "blended_pairwise" "away" "blended_pairwise" "blended" "away" "blended_pairwise" "blended")
LINESEARCHES=("agnostic" "agnostic" "adaptive" "adaptive" "adaptive" "adaptive" "adaptive" "adaptive")

N_COMBINATIONS=${#NORMS[@]}

# Combined CSV and summary
COMBINED_CSV="$BASE_DIR/run_selected/comparison_selected.csv"
BENCHMARK_SUMMARY="$BASE_DIR/run_selected/benchmark_summary_selected.txt"

# Accumulators for benchmark summary
declare -a SUMMARY_NAMES=()
declare -a SUMMARY_SUCCESS=()
declare -a SUMMARY_FAILED=()
declare -a SUMMARY_TOTAL=()
declare -A INSTANCE_BINVARS=()

grand_success=0
grand_failed=0
grand_total=0

grand_scip_timelimit=0
grand_global_timelimit=0
grand_restart_limit=0
grand_fw_infeasible=0
grand_other=0

grand_found_binary=0
grand_found_ginteger=0
grand_failed_binary=0
grand_failed_ginteger=0
grand_rr_found=0

# Clean up previous outputs before re-running
for i in "${!NORMS[@]}"; do
    rm -rf "$BASE_DIR/run_selected/${NORMS[$i]}_${VARIANTS[$i]}_${LINESEARCHES[$i]}_presolve_${PRESOLVE}/result"
done

CSV_HEADER="ID,Instance,BinaryVars,IntegerVars,SolutionFound,TotalTime,FPIterations,FWIterations,Restarts,Objective,Gap,FailureReason,FailureType,ProjectionNorm,FWVariant,LineSearch,RoundingThreshold,Presolve"
echo "$CSV_HEADER" > "$COMBINED_CSV"

for i in "${!NORMS[@]}"; do
    NORM="${NORMS[$i]}"
    VARIANT="${VARIANTS[$i]}"
    LS="${LINESEARCHES[$i]}"
    NAME="${NORM}_${VARIANT}_${LS}_presolve_${PRESOLVE}"

    OUTPUT_DIR="$BASE_DIR/run_selected/$NAME/output"
    RESULT_DIR="$BASE_DIR/run_selected/$NAME/result"
    mkdir -p "$RESULT_DIR"

    SUMMARY_FILE="${RESULT_DIR}/solution_summary.txt"
    DETAILED_FILE="${RESULT_DIR}/detailed_results.txt"
    FOUND_FILE="${RESULT_DIR}/solutions_found.txt"
    FAILED_FILE="${RESULT_DIR}/failed_runs.txt"
    CSV_FILE="${RESULT_DIR}/results.csv"

    # Initialize files
    echo "FPFW Solution Analysis [$NAME] - $(date)" > "$SUMMARY_FILE"
    echo "==========================================================" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"

    rounding_threshold=$(grep "Rounding thresh:" "$(ls "$OUTPUT_DIR"/job_*_*.out 2>/dev/null | head -1)" 2>/dev/null | head -1 | awk '{print $NF}')

    echo "DETAILED RESULTS" > "$DETAILED_FILE"
    echo "==========================================================" >> "$DETAILED_FILE"
    echo "  Projection norm:    $NORM" >> "$DETAILED_FILE"
    echo "  FW variant:         $VARIANT" >> "$DETAILED_FILE"
    echo "  Line search:        $LS" >> "$DETAILED_FILE"
    echo "  Rounding threshold: $rounding_threshold" >> "$DETAILED_FILE"
    echo "  Presolve:           $PRESOLVE" >> "$DETAILED_FILE"
    echo "==========================================================" >> "$DETAILED_FILE"
    echo "" >> "$DETAILED_FILE"

    echo "Test cases with solutions found:" > "$FOUND_FILE"
    echo "====================================================================================================================" >> "$FOUND_FILE"
    printf "%-3s %-25s %-8s %-12s %-10s %-10s %-10s %-12s %-8s\n" "ID" "Instance" "Bin Vars" "Time (s)" "FP Iters" "FW Iters" "Restarts" "Objective" "Gap (%)" >> "$FOUND_FILE"
    echo "====================================================================================================================" >> "$FOUND_FILE"

    echo "Failed/Interrupted runs:" > "$FAILED_FILE"
    echo "==========================================================================================================" >> "$FAILED_FILE"
    printf "%-3s %-25s %-8s %-15s %-50s\n" "ID" "Instance" "Bin Vars" "Runtime (s)" "Failure Reason" >> "$FAILED_FILE"
    echo "==========================================================================================================" >> "$FAILED_FILE"

    # CSV header
    echo "$CSV_HEADER" > "$CSV_FILE"

    # Counters
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

    # Statistics arrays
    declare -a found_times=()
    declare -a found_fw_iters=()
    declare -a found_fp_iters=()
    declare -a found_restarts=()

    # Function to extract job and task ID from filename
    extract_ids() {
        local filename=$(basename "$1")
        local job_id=$(echo "$filename" | sed -n 's/job_\([0-9]*\)_\([0-9]*\)\.\(out\|err\)/\1/p')
        local task_id=$(echo "$filename" | sed -n 's/job_\([0-9]*\)_\([0-9]*\)\.\(out\|err\)/\2/p')
        echo "$job_id $task_id"
    }

    # Process each output file
    for output_file in "$OUTPUT_DIR"/job_*_*.out; do
        if [ ! -f "$output_file" ]; then
            continue
        fi

        total_count=$((total_count + 1))

        read job_id task_id <<< $(extract_ids "$output_file")

        instance=$(grep "Running instance:" "$output_file" | head -1 | awk '{print $NF}')
        instance_name=$(basename "$instance" | sed 's/.mps.gz//' | sed 's/.mps//')

        projection_norm=$(grep "Projection norm:" "$output_file" | head -1 | awk '{print $NF}')
        if [ -z "$projection_norm" ]; then projection_norm="$NORM"; fi
        fw_variant=$(grep "FW variant:" "$output_file" | head -1 | awk '{print $NF}')
        line_search=$(grep "Line search:" "$output_file" | head -1 | awk '{print $NF}')
        rounding_threshold=$(grep "Rounding thresh:" "$output_file" | head -1 | awk '{print $NF}')

        # Parse exit reason
        exit_reason_str=$(grep "Exit reason:" "$output_file" | tail -1 | sed 's/Exit reason:[[:space:]]*//')

        scip_timelimit=""
        global_timelimit=""
        restart_limit=""
        fw_infeasible=""
        scip_solved=""
        rr_found=""

        case "$exit_reason_str" in
            *"SCIP time limit"*)               scip_timelimit="yes" ;;
            *"global time limit"*)             global_timelimit="yes" ;;
            *"cycled"*)                        restart_limit="yes" ;;
            *"outside the feasible polytope"*) fw_infeasible="yes" ;;
            *"solved by SCIP"*)                scip_solved="yes" ;;
            *"randomized rounding"*)           rr_found="yes" ;;
        esac

        # Extract performance metrics
        binary_vars=$(grep "Binary variables:" "$output_file" | tail -1 | awk '{print $NF}')
        integer_vars=$(grep "General integer variables:" "$output_file" | tail -1 | awk '{print $NF}')
        INSTANCE_BINVARS["$instance_name"]="$binary_vars"
        total_time=$(grep "Total time:" "$output_file" | tail -1 | awk '{print $3}')
        fw_time=$(grep "FW time:" "$output_file" | tail -1 | awk '{print $3}')
        fp_iterations=$(grep "FP iterations:" "$output_file" | tail -1 | awk '{print $NF}')
        fw_iterations=$(grep "FW iterations:" "$output_file" | tail -1 | awk '{print $NF}')
        restarts=$(grep "Restarts:" "$output_file" | tail -1 | awk '{print $NF}')

        solution_found=$(grep "Solution found:" "$output_file" | tail -1 | awk '{print $NF}')
        if [ -z "$solution_found" ]; then solution_found="false"; fi

        if [ "$solution_found" = "true" ]; then
            objective=$(grep "Objective:" "$output_file" | tail -1 | awk '{print $NF}')
            gap_line=$(grep "^Gap:" "$output_file" | tail -1)
            if echo "$gap_line" | grep -q "Infinite"; then
                gap="Infinite"
            else
                gap=$(echo "$gap_line" | awk '{print $(NF-1)}' | tr -d '%')
            fi
        else
            objective="N/A"
            gap="N/A"
        fi

        # Determine failure reason
        failure_reason="None"
        failure_type="None"
        this_was_found=0

        if [ -n "$scip_timelimit" ]; then
            failure_reason="SCIP_TIME_LIMIT (${DEF_SCIP_TIME_LIMIT}s solver limit)"
            failure_type="SCIP_TIMELIMIT"
            failed_count=$((failed_count + 1))
            scip_timelimit_count=$((scip_timelimit_count + 1))
            printf "%-3s %-25s %-8s %-15s %-50s\n" "$task_id" "$instance_name" "$binary_vars" "${total_time:-N/A}" "SCIP time limit (${DEF_SCIP_TIME_LIMIT}s)" >> "$FAILED_FILE"

        elif [ -n "$global_timelimit" ]; then
            failure_reason="GLOBAL_TIME_LIMIT (${DEF_GLOBAL_TIME_LIMIT}s per iteration)"
            failure_type="GLOBAL_TIMELIMIT"
            failed_count=$((failed_count + 1))
            global_timelimit_count=$((global_timelimit_count + 1))
            printf "%-3s %-25s %-8s %-15s %-50s\n" "$task_id" "$instance_name" "$binary_vars" "${total_time:-N/A}" "Global time limit (${DEF_GLOBAL_TIME_LIMIT}s)" >> "$FAILED_FILE"

        elif [ "$solution_found" = "true" ]; then
            found_count=$((found_count + 1))
            this_was_found=1
            [ -n "$rr_found" ] && rr_found_count=$((rr_found_count + 1))
            found_times+=("$total_time")
            found_fw_iters+=("$fw_iterations")
            found_fp_iters+=("$fp_iterations")
            found_restarts+=("$restarts")
            printf "%-3s %-25s %-8s %-12s %-10s %-10s %-10s %-12s %-8s\n" \
                "$task_id" "$instance_name" "$binary_vars" "$total_time" "$fp_iterations" "$fw_iterations" "$restarts" "$objective" "$gap" >> "$FOUND_FILE"

        elif [ -n "$fw_infeasible" ]; then
            failure_reason="FW_INFEASIBLE (FW returned infeasible solution)"
            failure_type="FW_INFEASIBLE"
            failed_count=$((failed_count + 1))
            fw_infeasible_count=$((fw_infeasible_count + 1))
            printf "%-3s %-25s %-8s %-15s %-50s\n" "$task_id" "$instance_name" "$binary_vars" "${total_time:-N/A}" "FW returns infeasible solution" >> "$FAILED_FILE"

        elif [ -n "$restart_limit" ]; then
            failure_reason="RESTART_LIMIT (max restarts reached)"
            failure_type="RESTART_LIMIT"
            failed_count=$((failed_count + 1))
            restart_limit_count=$((restart_limit_count + 1))
            printf "%-3s %-25s %-8s %-15s %-50s\n" "$task_id" "$instance_name" "$binary_vars" "${total_time:-N/A}" "Maximum restarts reached" >> "$FAILED_FILE"

        else
            failure_reason="UNKNOWN_ERROR"
            failure_type="UNKNOWN"
            failed_count=$((failed_count + 1))
            other_failure_count=$((other_failure_count + 1))
            printf "%-3s %-25s %-8s %-15s %-50s\n" "$task_id" "$instance_name" "$binary_vars" "${total_time:-N/A}" "Unknown error" >> "$FAILED_FILE"
        fi

        # Classify by instance type (binary-only vs general integer)
        if [ "${integer_vars:-0}" -gt 0 ] 2>/dev/null; then
            [ "$this_was_found" = "1" ] && found_ginteger=$((found_ginteger + 1)) || failed_ginteger=$((failed_ginteger + 1))
        else
            [ "$this_was_found" = "1" ] && found_binary=$((found_binary + 1)) || failed_binary=$((failed_binary + 1))
        fi

        # Write to CSVs
        echo "${task_id},${instance_name},${binary_vars},${integer_vars},${solution_found},${total_time},${fp_iterations},${fw_iterations},${restarts},${objective},${gap},${failure_reason},${failure_type},${projection_norm},${fw_variant},${line_search},${rounding_threshold},${PRESOLVE}" >> "$CSV_FILE"
        echo "${task_id},${instance_name},${binary_vars},${integer_vars},${solution_found},${total_time},${fp_iterations},${fw_iterations},${restarts},${objective},${gap},${failure_reason},${failure_type},${projection_norm},${fw_variant},${line_search},${rounding_threshold},${PRESOLVE}" >> "$COMBINED_CSV"

        # Write to detailed file
        echo "Instance: ${instance_name} (ID:${task_id})" >> "$DETAILED_FILE"
        echo "  Binary variables:  ${binary_vars}" >> "$DETAILED_FILE"
        echo "  Integer variables: ${integer_vars}" >> "$DETAILED_FILE"
        echo "  Solution found:    ${solution_found}" >> "$DETAILED_FILE"
        echo "  Total time:        ${total_time}" >> "$DETAILED_FILE"
        echo "  FP iterations:     ${fp_iterations}" >> "$DETAILED_FILE"
        echo "  FW iterations:     ${fw_iterations}" >> "$DETAILED_FILE"
        echo "  FW time:           ${fw_time}" >> "$DETAILED_FILE"
        echo "  Restarts:          ${restarts}" >> "$DETAILED_FILE"

        if [ "$solution_found" = "true" ]; then
            echo "  Objective:         ${objective}" >> "$DETAILED_FILE"
            if [ "$gap" = "Infinite" ]; then
                echo "  Gap:               Infinite" >> "$DETAILED_FILE"
            else
                echo "  Gap:               ${gap}%" >> "$DETAILED_FILE"
            fi
        else
            echo "  Failure reason:    ${failure_reason}" >> "$DETAILED_FILE"
            echo "  Failure type:      ${failure_type}" >> "$DETAILED_FILE"
        fi
        echo "" >> "$DETAILED_FILE"
    done

    # Calculate statistics for successful runs
    if [ ${#found_times[@]} -gt 0 ]; then
        sum_time=0
        for t in "${found_times[@]}"; do
            sum_time=$(awk -v s="$sum_time" -v v="$t" 'BEGIN {print s + v}')
        done
        avg_time=$(awk -v s="$sum_time" -v n="${#found_times[@]}" 'BEGIN {printf "%.2f", s / n}')

        sum_fw=0
        for i in "${found_fw_iters[@]}"; do sum_fw=$((sum_fw + i)); done
        avg_fw=$((sum_fw / ${#found_fw_iters[@]}))

        sum_fp=0
        for i in "${found_fp_iters[@]}"; do sum_fp=$((sum_fp + i)); done
        avg_fp=$((sum_fp / ${#found_fp_iters[@]}))

        sum_restarts=0
        for i in "${found_restarts[@]}"; do sum_restarts=$((sum_restarts + i)); done
        avg_restarts=$(awk -v s="$sum_restarts" -v n="${#found_restarts[@]}" 'BEGIN {printf "%.2f", s / n}')
    fi

    # Write summary
    echo "OVERALL STATISTICS" >> "$SUMMARY_FILE"
    echo "==========================================================" >> "$SUMMARY_FILE"
    echo "Total test cases processed: $total_count" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    echo "Solutions found: $found_count ($(awk -v f="$found_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)" >> "$SUMMARY_FILE"
    echo "  - Via randomized rounding: $rr_found_count" >> "$SUMMARY_FILE"
    echo "  - Via FPFW main loop:      $((found_count - rr_found_count))" >> "$SUMMARY_FILE"
    echo "Failed/Interrupted: $failed_count ($(awk -v f="$failed_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"

    echo "BREAKDOWN BY INSTANCE TYPE" >> "$SUMMARY_FILE"
    echo "==========================================================" >> "$SUMMARY_FILE"
    echo "Binary-only (20):      found=$found_binary  failed=$failed_binary" >> "$SUMMARY_FILE"
    echo "General integer (5):   found=$found_ginteger  failed=$failed_ginteger" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"

    echo "FAILURE BREAKDOWN" >> "$SUMMARY_FILE"
    echo "==========================================================" >> "$SUMMARY_FILE"
    echo "Total failed cases: $failed_count" >> "$SUMMARY_FILE"
    echo "  - Due to SCIP time limit:  $scip_timelimit_count" >> "$SUMMARY_FILE"
    echo "  - Due to global time limit: $global_timelimit_count" >> "$SUMMARY_FILE"
    echo "  - Due to restart limit:    $restart_limit_count" >> "$SUMMARY_FILE"
    echo "  - Due to FW infeasibility: $fw_infeasible_count" >> "$SUMMARY_FILE"
    echo "  - Due to other reasons:    $other_failure_count" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"

    if [ ${#found_times[@]} -gt 0 ]; then
        echo "STATISTICS FOR SUCCESSFUL RUNS" >> "$SUMMARY_FILE"
        echo "==========================================================" >> "$SUMMARY_FILE"
        echo "Average time:          ${avg_time}s" >> "$SUMMARY_FILE"
        echo "Average FP iterations: $avg_fp" >> "$SUMMARY_FILE"
        echo "Average FW iterations: $avg_fw" >> "$SUMMARY_FILE"
        echo "Average restarts:      $avg_restarts" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
    fi

    echo "OUTPUT FILES" >> "$SUMMARY_FILE"
    echo "==========================================================" >> "$SUMMARY_FILE"
    echo "  - Summary:  $SUMMARY_FILE" >> "$SUMMARY_FILE"
    echo "  - Detailed: $DETAILED_FILE" >> "$SUMMARY_FILE"
    echo "  - Found:    $FOUND_FILE" >> "$SUMMARY_FILE"
    echo "  - Failed:   $FAILED_FILE" >> "$SUMMARY_FILE"
    echo "  - CSV:      $CSV_FILE" >> "$SUMMARY_FILE"

    echo ""
    echo "=========================================================="
    cat "$SUMMARY_FILE"
    echo "=========================================================="
    echo ""
    echo "Results saved to $RESULT_DIR"

    # Accumulate for benchmark summary
    SUMMARY_NAMES+=("$NAME")
    SUMMARY_SUCCESS+=("$found_count")
    SUMMARY_FAILED+=("$failed_count")
    SUMMARY_TOTAL+=("$total_count")
    grand_success=$((grand_success + found_count))
    grand_failed=$((grand_failed + failed_count))
    grand_total=$((grand_total + total_count))
    grand_scip_timelimit=$((grand_scip_timelimit + scip_timelimit_count))
    grand_global_timelimit=$((grand_global_timelimit + global_timelimit_count))
    grand_restart_limit=$((grand_restart_limit + restart_limit_count))
    grand_fw_infeasible=$((grand_fw_infeasible + fw_infeasible_count))
    grand_other=$((grand_other + other_failure_count))
    grand_found_binary=$((grand_found_binary + found_binary))
    grand_found_ginteger=$((grand_found_ginteger + found_ginteger))
    grand_failed_binary=$((grand_failed_binary + failed_binary))
    grand_failed_ginteger=$((grand_failed_ginteger + failed_ginteger))
    grand_rr_found=$((grand_rr_found + rr_found_count))
done

# ── Generate benchmark_summary_selected.txt ───────────────────────────────────
{
echo "FPFW Benchmark Summary (Selected Combinations)"
echo "Generated: $(date)"
echo "Instances: $(( grand_total / N_COMBINATIONS )) per combination | Total combinations: $N_COMBINATIONS"
echo ""
printf "%-65s %7s %7s %7s %7s\n" "Combination (norm_variant_linesearch_presolve)" "Success" "Failed" "Total" "Rate(%)"
echo "------------------------------------------------------------------------------------------------"
for j in "${!SUMMARY_NAMES[@]}"; do
    rate=$(awk -v s="${SUMMARY_SUCCESS[$j]}" -v t="${SUMMARY_TOTAL[$j]}" \
        'BEGIN {if(t>0) printf "%.1f", (s/t)*100; else print "0.0"}')
    printf "%-65s %7d %7d %7d %7s\n" \
        "${SUMMARY_NAMES[$j]}" "${SUMMARY_SUCCESS[$j]}" "${SUMMARY_FAILED[$j]}" "${SUMMARY_TOTAL[$j]}" "$rate"
done
overall_rate=$(awk -v s="$grand_success" -v t="$grand_total" \
    'BEGIN {if(t>0) printf "%.1f", (s/t)*100; else print "0.0"}')
echo "------------------------------------------------------------------------------------------------"
printf "%-65s %7d %7d %7d %7s\n" "TOTAL" "$grand_success" "$grand_failed" "$grand_total" "$overall_rate"
echo ""
echo "SOLUTION METHOD BREAKDOWN (across all combinations)"
echo "-----------------------------------------------------"
echo "  Via randomized rounding: $grand_rr_found"
echo "  Via FPFW main loop:      $((grand_success - grand_rr_found))"
echo ""
echo "BREAKDOWN BY INSTANCE TYPE (across all combinations)"
echo "-----------------------------------------------------"
echo "  Binary-only (20):      found=$grand_found_binary  failed=$grand_failed_binary"
echo "  General integer (5):   found=$grand_found_ginteger  failed=$grand_failed_ginteger"
echo ""
echo "FAILURE BREAKDOWN (across all combinations)"
echo "--------------------------------------------"
echo "  Global time limit: $grand_global_timelimit"
echo "  SCIP time limit:   $grand_scip_timelimit"
echo "  Restart limit:     $grand_restart_limit"
echo "  FW infeasible:     $grand_fw_infeasible"
echo "  Unknown:           $grand_other"
echo ""
echo "INSTANCE BINARY VARIABLE COUNTS"
echo "------------------------------------"
printf "%-40s %s\n" "Instance" "Bin Vars"
echo "------------------------------------"
for inst in $(echo "${!INSTANCE_BINVARS[@]}" | tr ' ' '\n' | sort); do
    printf "%-40s %s\n" "$inst" "${INSTANCE_BINVARS[$inst]}"
done
echo "------------------------------------"
} | tee "$BENCHMARK_SUMMARY"

echo ""
echo "All $N_COMBINATIONS combinations checked."
echo "Combined CSV:       $COMBINED_CSV"
echo "Benchmark summary:  $BENCHMARK_SUMMARY"
