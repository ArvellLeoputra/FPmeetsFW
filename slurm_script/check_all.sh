#!/bin/bash

# Check results for all 28 FP-FW experiment combinations
# Generates individual result files per combination + one combined comparison CSV

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"

NORMS=("manhattan" "euclidean" "abssmooth" "euclidean" "euclidean" "abssmooth" "abssmooth")
LINESEARCHES=("agnostic" "agnostic" "agnostic" "adaptive" "secant" "adaptive" "secant")
VARIANTS=("vanilla" "away" "blended_pairwise" "blended")

# Combined CSV across all 28 runs
COMBINED_CSV="$BASE_DIR/run/comparison.csv"
echo "TaskID,Instance,Norm,Variant,LineSearch,Status,TotalTime,FPIterations,FWIterations,RestartsCycles,RestartsFixed,SolutionFound,Objective,Gap,FailureReason" > "$COMBINED_CSV"

for i in "${!NORMS[@]}"; do
    NORM="${NORMS[$i]}"
    LS="${LINESEARCHES[$i]}"

    for VARIANT in "${VARIANTS[@]}"; do
        NAME="${NORM}_${VARIANT}_${LS}"
        OUTPUT_DIR="$BASE_DIR/run/$NAME/output"
        RESULT_DIR="$BASE_DIR/run/$NAME/result"
        mkdir -p "$RESULT_DIR"

        SUMMARY_FILE="$RESULT_DIR/solution_summary.txt"
        DETAILED_FILE="$RESULT_DIR/detailed_results.txt"
        FOUND_FILE="$RESULT_DIR/solutions_found.txt"
        NOT_FOUND_FILE="$RESULT_DIR/solutions_not_found.txt"
        FAILED_FILE="$RESULT_DIR/failed_runs.txt"
        CSV_FILE="$RESULT_DIR/results.csv"

        # Initialize files
        echo "FPFW Solution Analysis [$NAME] - $(date)" > "$SUMMARY_FILE"
        echo "==========================================================" >> "$SUMMARY_FILE"

        echo "DETAILED RESULTS" > "$DETAILED_FILE"
        echo "==========================================================" >> "$DETAILED_FILE"

        echo "Test cases with solutions found:" > "$FOUND_FILE"
        printf "%-8s %-40s %-12s %-10s %-10s %-8s %-12s\n" "Task ID" "Instance" "Time (s)" "FP Iters" "FW Iters" "Gap (%)" "Objective" >> "$FOUND_FILE"
        echo "==========================================================" >> "$FOUND_FILE"

        echo "Test cases without solutions:" > "$NOT_FOUND_FILE"
        printf "%-8s %-40s %-12s %-10s %-10s %-40s\n" "Task ID" "Instance" "Time (s)" "FP Iters" "FW Iters" "Reason" >> "$NOT_FOUND_FILE"
        echo "==========================================================" >> "$NOT_FOUND_FILE"

        echo "Failed/Interrupted runs:" > "$FAILED_FILE"
        printf "%-8s %-40s %-15s %-50s\n" "Task ID" "Instance" "Runtime (s)" "Failure Reason" >> "$FAILED_FILE"
        echo "==========================================================" >> "$FAILED_FILE"

        echo "TaskID,Instance,Status,TotalTime,FPIterations,FWIterations,RestartsCycles,RestartsFixed,SolutionFound,Objective,Gap,FailureReason" > "$CSV_FILE"

        total_count=0
        found_count=0
        failed_count=0
        scip_timelimit_count=0
        fw_timelimit_count=0
        restart_limit_count=0
        fw_infeasible_count=0
        other_failure_count=0

        declare -a found_times=()
        declare -a found_fw_iters=()
        declare -a found_fp_iters=()

        for output_file in "$OUTPUT_DIR"/job_*_*.out; do
            [ -f "$output_file" ] || continue
            total_count=$((total_count + 1))

            filename=$(basename "$output_file")
            task_id=$(echo "$filename" | sed -n 's/job_[0-9]*_\([0-9]*\)\.out/\1/p')

            instance=$(grep "Running instance:" "$output_file" | head -1 | awk '{print $NF}')
            instance_name=$(basename "$instance" | sed 's/.mps.gz//' | sed 's/.mps//')

            total_time=$(grep "Total time:" "$output_file" | tail -1 | awk '{print $3}')
            fp_iterations=$(grep "FP iterations:" "$output_file" | tail -1 | awk '{print $NF}')
            fw_iterations=$(grep "FW iterations:" "$output_file" | tail -1 | awk '{print $NF}')
            restarts_cycles=$(grep "Restarts (cycles):" "$output_file" | tail -1 | awk '{print $NF}')
            restarts_fixed=$(grep "Restarts (fixed):" "$output_file" | tail -1 | awk '{print $NF}')
            solution_found=$(grep "Solution found:" "$output_file" | tail -1 | awk '{print $NF}')
            [ -z "$solution_found" ] && solution_found="false"

            scip_timelimit=$(grep -q "time limit reached" "$output_file" && echo "yes" || echo "")
            fw_timelimit=$(grep -q "FW time limit reached" "$output_file" && echo "yes" || echo "")
            cycle_restart_limit=$(grep -q "Maximum cycle restarts.*reached" "$output_file" && echo "yes" || echo "")
            fixed_restart_limit=$(grep -q "Maximum fixed-point restarts.*reached" "$output_file" && echo "yes" || echo "")
            restart_limit=""
            [ -n "$cycle_restart_limit" ] || [ -n "$fixed_restart_limit" ] && restart_limit="yes"
            fw_infeasible=$(grep -q "FrankWolfe return infeasible solution" "$output_file" && echo "yes" || echo "")

            objective="N/A"
            gap="N/A"
            if [ "$solution_found" = "true" ]; then
                objective=$(grep "Objective:" "$output_file" | tail -1 | awk '{print $NF}')
                gap=$(grep "Gap" "$output_file" | grep "%" | tail -1 | awk '{print $(NF-1)}' | tr -d '%')
            fi

            status="UNKNOWN"
            failure_reason="None"

            if [ -n "$scip_timelimit" ]; then
                status="FAILED"; failure_reason="SCIP_TIME_LIMIT"
                failed_count=$((failed_count + 1)); scip_timelimit_count=$((scip_timelimit_count + 1))
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "${total_time:-N/A}" "SCIP time limit" >> "$FAILED_FILE"
            elif [ -n "$fw_timelimit" ]; then
                status="FAILED"; failure_reason="FW_TIME_LIMIT"
                failed_count=$((failed_count + 1)); fw_timelimit_count=$((fw_timelimit_count + 1))
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "${total_time:-N/A}" "FW time limit" >> "$FAILED_FILE"
            elif [ "$solution_found" = "true" ]; then
                status="SUCCESS"
                found_count=$((found_count + 1))
                found_times+=("$total_time")
                found_fw_iters+=("$fw_iterations")
                found_fp_iters+=("$fp_iterations")
                printf "%-8s %-40s %-12s %-10s %-10s %-8s %-12s\n" \
                    "$task_id" "$instance_name" "$total_time" "$fp_iterations" "$fw_iterations" "$gap" "$objective" >> "$FOUND_FILE"
            elif [ -n "$fw_infeasible" ]; then
                status="FAILED"; failure_reason="FW_INFEASIBLE"
                failed_count=$((failed_count + 1)); fw_infeasible_count=$((fw_infeasible_count + 1))
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "${total_time:-N/A}" "FW infeasible" >> "$FAILED_FILE"
            elif [ -n "$restart_limit" ]; then
                status="FAILED"; failure_reason="RESTART_LIMIT"
                failed_count=$((failed_count + 1)); restart_limit_count=$((restart_limit_count + 1))
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "${total_time:-N/A}" "Restart limit reached" >> "$FAILED_FILE"
            else
                status="FAILED"; failure_reason="UNKNOWN"
                failed_count=$((failed_count + 1)); other_failure_count=$((other_failure_count + 1))
                printf "%-8s %-40s %-15s %-50s\n" "$task_id" "$instance_name" "${total_time:-N/A}" "Unknown" >> "$FAILED_FILE"
            fi

            echo "${task_id},${instance_name},${status},${total_time},${fp_iterations},${fw_iterations},${restarts_cycles},${restarts_fixed},${solution_found},${objective},${gap},${failure_reason}" >> "$CSV_FILE"
            echo "${task_id},${instance_name},${NORM},${VARIANT},${LS},${status},${total_time},${fp_iterations},${fw_iterations},${restarts_cycles},${restarts_fixed},${solution_found},${objective},${gap},${failure_reason}" >> "$COMBINED_CSV"

            echo "Task ${task_id}: ${instance_name}" >> "$DETAILED_FILE"
            echo "  Status: ${status}" >> "$DETAILED_FILE"
            echo "  Total time: ${total_time}s" >> "$DETAILED_FILE"
            echo "  FP iterations: ${fp_iterations}" >> "$DETAILED_FILE"
            echo "  FW iterations: ${fw_iterations}" >> "$DETAILED_FILE"
            echo "  Restarts (cycles): ${restarts_cycles}" >> "$DETAILED_FILE"
            echo "  Restarts (fixed): ${restarts_fixed}" >> "$DETAILED_FILE"
            echo "  Solution found: ${solution_found}" >> "$DETAILED_FILE"
            [ "$solution_found" = "true" ] && echo "  Objective: ${objective}" >> "$DETAILED_FILE"
            [ "$solution_found" = "true" ] && echo "  Gap: ${gap}%" >> "$DETAILED_FILE"
            [ "$status" = "FAILED" ] && echo "  Failure reason: ${failure_reason}" >> "$DETAILED_FILE"
            echo "" >> "$DETAILED_FILE"
        done

        {
            echo ""
            echo "OVERALL STATISTICS"
            echo "=========================================================="
            echo "Total processed: $total_count"
            echo "Solutions found: $found_count"
            echo "Failed:          $failed_count"
            echo ""
            echo "FAILURE BREAKDOWN"
            echo "  SCIP time limit:  $scip_timelimit_count"
            echo "  FW time limit:    $fw_timelimit_count"
            echo "  Restart limit:    $restart_limit_count"
            echo "  FW infeasible:    $fw_infeasible_count"
            echo "  Other:            $other_failure_count"
        } >> "$SUMMARY_FILE"

        echo "Done: $NAME ($found_count/$total_count solved)"
    done
done

echo ""
echo "All 28 combinations checked."
echo "Combined CSV: $COMBINED_CSV"
