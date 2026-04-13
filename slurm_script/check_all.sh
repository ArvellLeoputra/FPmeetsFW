#!/bin/bash

# Check results for all 28 FP-FW experiment combinations
# Also generates one combined comparison CSV at the end

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"

NORMS=("manhattan" "euclidean" "abssmooth" "euclidean" "euclidean" "abssmooth" "abssmooth")
LINESEARCHES=("agnostic" "agnostic" "agnostic" "adaptive" "secant" "adaptive" "secant")
VARIANTS=("vanilla" "away" "blended_pairwise" "blended")

# Combined CSV across all 28 runs
COMBINED_CSV="$BASE_DIR/run/comparison.csv"

# Clean up previous check_all.sh outputs before re-running
rm -f "$COMBINED_CSV"
for i in "${!NORMS[@]}"; do
    for VARIANT in "${VARIANTS[@]}"; do
        rm -rf "$BASE_DIR/run/${NORMS[$i]}_${VARIANT}_${LINESEARCHES[$i]}/result"
    done
done

echo "TaskID,Instance,Status,TotalTime,FPIterations,FWIterations,RestartsCycles,RestartsFixed,SolutionFound,Objective,Gap,FailureReason,FailureType,ProjectionNorm,FWVariant,LineSearch,RoundingThreshold" > "$COMBINED_CSV"

for i in "${!NORMS[@]}"; do
    NORM="${NORMS[$i]}"
    LS="${LINESEARCHES[$i]}"

    for VARIANT in "${VARIANTS[@]}"; do
        NAME="${NORM}_${VARIANT}_${LS}"

        OUTPUT_DIR="$BASE_DIR/run/$NAME/output"
        ERROR_DIR="$BASE_DIR/run/$NAME/error"
        RESULT_DIR="$BASE_DIR/run/$NAME/result"
        mkdir -p "$RESULT_DIR"

        SUMMARY_FILE="${RESULT_DIR}/solution_summary.txt"
        DETAILED_FILE="${RESULT_DIR}/detailed_results.txt"
        FOUND_FILE="${RESULT_DIR}/solutions_found.txt"
        NOT_FOUND_FILE="${RESULT_DIR}/solutions_not_found.txt"
        FAILED_FILE="${RESULT_DIR}/failed_runs.txt"
        CSV_FILE="${RESULT_DIR}/results.csv"

        # Initialize files
        echo "FPFW Solution Analysis [$NAME] - $(date)" > "$SUMMARY_FILE"
        echo "==========================================================" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"

        echo "DETAILED RESULTS" > "$DETAILED_FILE"
        echo "==========================================================" >> "$DETAILED_FILE"
        echo "" >> "$DETAILED_FILE"

        echo "Test cases with solutions found:" > "$FOUND_FILE"
        echo "==========================================================" >> "$FOUND_FILE"
        printf "%-8s %-40s %-12s %-10s %-10s %-8s %-12s\n" "Task ID" "Instance" "Time (s)" "FP Iters" "FW Iters" "Gap (%)" "Objective" >> "$FOUND_FILE"
        echo "==========================================================" >> "$FOUND_FILE"

        echo "Test cases without solutions:" > "$NOT_FOUND_FILE"
        echo "==========================================================" >> "$NOT_FOUND_FILE"
        printf "%-8s %-40s %-12s %-10s %-10s %-40s\n" "Task ID" "Instance" "Time (s)" "FP Iters" "FW Iters" "Reason" >> "$NOT_FOUND_FILE"
        echo "==========================================================" >> "$NOT_FOUND_FILE"

        echo "Failed/Interrupted runs:" > "$FAILED_FILE"
        echo "==========================================================" >> "$FAILED_FILE"
        printf "%-8s %-40s %-15s %-50s\n" "Task ID" "Instance" "Runtime (s)" "Failure Reason" >> "$FAILED_FILE"
        echo "==========================================================" >> "$FAILED_FILE"

        # CSV header
        echo "TaskID,Instance,Status,TotalTime,FPIterations,FWIterations,RestartsCycles,RestartsFixed,SolutionFound,Objective,Gap,FailureReason,FailureType,ProjectionNorm,FWVariant,LineSearch,RoundingThreshold" > "$CSV_FILE"

        # Counters
        total_count=0
        found_count=0
        not_found_count=0
        failed_count=0
        slurm_timelimit_count=0
        scip_timelimit_count=0
        fw_timelimit_count=0
        restart_limit_count=0
        fw_infeasible_count=0
        memory_failure_count=0
        numerical_error_count=0
        other_failure_count=0
        fp_iterationlimit_count=0

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

            # Extract job and task IDs
            read job_id task_id <<< $(extract_ids "$output_file")

            # Extract basic information from output file
            instance=$(grep "Running instance:" "$output_file" | head -1 | awk '{print $NF}')
            instance_name=$(basename "$instance" | sed 's/.mps.gz//' | sed 's/.mps//')

            # Extract projection norm, FW variant, line search, rounding threshold
            projection_norm=$(grep "Projection norm:" "$output_file" | head -1 | awk '{print $NF}')
            if [ -z "$projection_norm" ]; then
                projection_norm="$NORM"
            fi
            fw_variant=$(grep "FW variant:" "$output_file" | head -1 | awk '{print $NF}')
            line_search=$(grep "Line search:" "$output_file" | head -1 | awk '{print $NF}')
            rounding_threshold=$(grep "Rounding threshold:" "$output_file" | head -1 | awk '{print $NF}')

            # No error files locally — skip SLURM/memory/numerical checks
            slurm_timelimit=""
            memory_failure=""
            numerical_error=""

            # Check for SCIP time limit
            scip_timelimit=""
            scip_status=$(grep "SCIP Status" "$output_file" | grep -i "time limit reached")
            if [ -n "$scip_status" ]; then
                scip_timelimit="yes"
            fi

            # Check for FW internal time limit
            fw_timelimit=$(grep -q "FW time limit reached" "$output_file" && echo "yes" || echo "")

            # Check for node limit
            node_limit=$(grep -q "node limit reached" "$output_file" && echo "yes" || echo "")

            # Check for restart limit (both cycle and fixed-point)
            restart_limit=$(grep -qE "Maximum (cycle|fixed-point) restarts.*reached, stopping" "$output_file" && echo "yes" || echo "")

            # Check for FW infeasibility
            fw_infeasible=$(grep -q "FrankWolfe return infeasible solution!" "$output_file" && echo "yes" || echo "")

            # Extract performance metrics
            total_time=$(grep "Total time:" "$output_file" | tail -1 | awk '{print $3}')
            fp_iterations=$(grep "FP iterations:" "$output_file" | tail -1 | awk '{print $NF}')
            fw_iterations=$(grep "FW iterations:" "$output_file" | tail -1 | awk '{print $NF}')
            restarts_cycles=$(grep "Restarts (cycles):" "$output_file" | tail -1 | awk '{print $NF}')
            restarts_fixed=$(grep "Restarts (fixed):" "$output_file" | tail -1 | awk '{print $NF}')

            # Check if solution was found
            solution_found=$(grep "Solution found:" "$output_file" | tail -1 | awk '{print $NF}')
            if [ -z "$solution_found" ]; then
                solution_found="false"
            fi

            # Extract solution details if found
            if [ "$solution_found" = "true" ]; then
                objective=$(grep "Objective:" "$output_file" | tail -1 | awk '{print $NF}')
                gap=$(grep "Gap" "$output_file" | grep "%" | tail -1 | awk '{print $(NF-1)}' | tr -d '%')
                scip_final_status=$(grep "Status:" "$output_file" | tail -1 | awk -F': ' '{print $2}')
            else
                objective="N/A"
                gap="N/A"
                scip_final_status="N/A"
            fi

            # Determine status and failure reason
            status="UNKNOWN"
            failure_reason="None"
            failure_type="None"

            # Priority 1: SCIP time limit
            if [ -n "$scip_timelimit" ]; then
                status="FAILED"
                failure_reason="SCIP_TIME_LIMIT (480s solver limit)"
                failure_type="SCIP_TIMELIMIT"
                failed_count=$((failed_count + 1))
                scip_timelimit_count=$((scip_timelimit_count + 1))
                runtime_display="${total_time:-N/A}"
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "$runtime_display" "SCIP time limit (480s)" >> "$FAILED_FILE"

            # Priority 2: FW time limit
            elif [ -n "$fw_timelimit" ]; then
                status="FAILED"
                failure_reason="FW_TIME_LIMIT (300s per iteration)"
                failure_type="FW_TIMELIMIT"
                failed_count=$((failed_count + 1))
                fw_timelimit_count=$((fw_timelimit_count + 1))
                runtime_display="${total_time:-N/A}"
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "$runtime_display" "FW time limit (300s)" >> "$FAILED_FILE"

            # Priority 3: Solution found
            elif [ "$solution_found" = "true" ]; then
                status="SUCCESS"
                found_count=$((found_count + 1))
                found_times+=("$total_time")
                found_fw_iters+=("$fw_iterations")
                found_fp_iters+=("$fp_iterations")
                found_restarts+=("$restarts_cycles")
                printf "%-8s %-40s %-12s %-10s %-10s %-10s %-8s %-12s\n" \
                    "$task_id" "$instance_name" "$total_time" "$fp_iterations" "$fw_iterations" "$restarts_cycles" "$gap" "$objective" >> "$FOUND_FILE"

            # Priority 4: FW returns infeasible
            elif [ -n "$fw_infeasible" ]; then
                status="FAILED"
                failure_reason="FW_INFEASIBLE (FW returned infeasible solution)"
                failure_type="FW_INFEASIBLE"
                failed_count=$((failed_count + 1))
                fw_infeasible_count=$((fw_infeasible_count + 1))
                runtime_display="${total_time:-N/A}"
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "$runtime_display" "FW returns infeasible solution" >> "$FAILED_FILE"

            # Priority 5: Restart limit
            elif [ -n "$restart_limit" ]; then
                status="FAILED"
                failure_reason="RESTART_LIMIT (max restarts reached)"
                failure_type="RESTART_LIMIT"
                failed_count=$((failed_count + 1))
                restart_limit_count=$((restart_limit_count + 1))
                runtime_display="${total_time:-N/A}"
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "$runtime_display" "Maximum restarts reached" >> "$FAILED_FILE"

            # Priority 6: Unknown
            else
                status="FAILED"
                failure_reason="UNKNOWN_ERROR"
                failure_type="UNKNOWN"
                failed_count=$((failed_count + 1))
                other_failure_count=$((other_failure_count + 1))
                runtime_display="${total_time:-N/A}"
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "$runtime_display" "Unknown error" >> "$FAILED_FILE"
            fi

            # Write to individual CSV
            echo "${task_id},${instance_name},${status},${total_time},${fp_iterations},${fw_iterations},${restarts_cycles},${restarts_fixed},${solution_found},${objective},${gap},${failure_reason},${failure_type},${projection_norm},${fw_variant},${line_search},${rounding_threshold}" >> "$CSV_FILE"

            # Write to combined CSV
            echo "${task_id},${instance_name},${status},${total_time},${fp_iterations},${fw_iterations},${restarts_cycles},${restarts_fixed},${solution_found},${objective},${gap},${failure_reason},${failure_type},${projection_norm},${fw_variant},${line_search},${rounding_threshold}" >> "$COMBINED_CSV"

            # Write to detailed file
            echo "Task ${task_id}: ${instance_name}" >> "$DETAILED_FILE"
            echo "  Status: ${status}" >> "$DETAILED_FILE"
            echo "  Projection norm: ${projection_norm}" >> "$DETAILED_FILE"
            echo "  FW variant: ${fw_variant}" >> "$DETAILED_FILE"
            echo "  Line search: ${line_search}" >> "$DETAILED_FILE"
            echo "  Rounding threshold: ${rounding_threshold}" >> "$DETAILED_FILE"

            if [ "$solution_found" = "true" ] || [ "$failure_type" = "RESTART_LIMIT" ] || [ "$failure_type" = "FW_TIMELIMIT" ]; then
                echo "  Total time: ${total_time}s" >> "$DETAILED_FILE"
                echo "  FP iterations: ${fp_iterations}" >> "$DETAILED_FILE"
                echo "  FW iterations: ${fw_iterations}" >> "$DETAILED_FILE"
                echo "  Restarts (cycles): ${restarts_cycles}" >> "$DETAILED_FILE"
                echo "  Restarts (fixed): ${restarts_fixed}" >> "$DETAILED_FILE"
            fi

            echo "  Solution found: ${solution_found}" >> "$DETAILED_FILE"

            if [ "$solution_found" = "true" ]; then
                echo "  Objective: ${objective}" >> "$DETAILED_FILE"
                echo "  Gap: ${gap}%" >> "$DETAILED_FILE"
            else
                echo "  Failure reason: ${failure_reason}" >> "$DETAILED_FILE"
                echo "  Failure type: ${failure_type}" >> "$DETAILED_FILE"
                if [ -n "$scip_timelimit" ]; then
                    echo "  Note: SCIP reached 480s time limit" >> "$DETAILED_FILE"
                fi
                if [ -n "$fw_infeasible" ]; then
                    echo "  Note: FrankWolfe returned infeasible solution" >> "$DETAILED_FILE"
                fi
                if [ -n "$restart_limit" ]; then
                    echo "  Note: Maximum restarts reached" >> "$DETAILED_FILE"
                fi
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
            for i in "${found_fw_iters[@]}"; do
                sum_fw=$((sum_fw + i))
            done
            avg_fw=$((sum_fw / ${#found_fw_iters[@]}))

            sum_fp=0
            for i in "${found_fp_iters[@]}"; do
                sum_fp=$((sum_fp + i))
            done
            avg_fp=$((sum_fp / ${#found_fp_iters[@]}))

            sum_restarts=0
            for i in "${found_restarts[@]}"; do
                sum_restarts=$((sum_restarts + i))
            done
            avg_restarts=$(awk -v s="$sum_restarts" -v n="${#found_restarts[@]}" 'BEGIN {printf "%.2f", s / n}')
        fi

        # Write summary
        echo "OVERALL STATISTICS" >> "$SUMMARY_FILE"
        echo "==========================================================" >> "$SUMMARY_FILE"
        echo "Total test cases processed: $total_count" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"
        echo "Solutions found: $found_count ($(awk -v f="$found_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)" >> "$SUMMARY_FILE"
        echo "Failed/Interrupted: $failed_count ($(awk -v f="$failed_count" -v t="$total_count" 'BEGIN {if(t>0) printf "%.2f", (f/t)*100; else print "0.00"}')%)" >> "$SUMMARY_FILE"
        echo "" >> "$SUMMARY_FILE"

        echo "FAILURE BREAKDOWN" >> "$SUMMARY_FILE"
        echo "==========================================================" >> "$SUMMARY_FILE"
        echo "Total failed cases: $failed_count" >> "$SUMMARY_FILE"
        echo "  - Due to SCIP time limit:  $scip_timelimit_count" >> "$SUMMARY_FILE"
        echo "  - Due to FW time limit:    $fw_timelimit_count" >> "$SUMMARY_FILE"
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

        # Print summary to console
        echo ""
        echo "=========================================================="
        cat "$SUMMARY_FILE"
        echo "=========================================================="
        echo ""
        echo "Results saved to $RESULT_DIR"
    done
done

echo ""
echo "All 28 combinations checked."
echo "Combined CSV: $COMBINED_CSV"
