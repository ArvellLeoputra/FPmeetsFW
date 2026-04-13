#!/bin/bash

# Run test1-3 with DEBUG_VERBOSE=true, then restore

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"
TESTCASE_DIR="$BASE_DIR/testcase"
DEP_FILE="$BASE_DIR/dependencies.jl"

# Flip DEBUG_VERBOSE to true
sed -i 's/const DEBUG_VERBOSE = false/const DEBUG_VERBOSE = true/' "$DEP_FILE"
echo "DEBUG_VERBOSE set to true"

for i in 1 2 3; do
    INSTANCE="$TESTCASE_DIR/test${i}.mps"
    OUTPUT="$BASE_DIR/testcase/test${i}_output.txt"
    echo "Running test${i}..."
    julia --project "$BASE_DIR/run_test.jl" "$INSTANCE" manhattan 0.5 vanilla agnostic > "$OUTPUT" 2>&1
    echo "Done. Output: $OUTPUT"
done

# Restore DEBUG_VERBOSE to false
sed -i 's/const DEBUG_VERBOSE = true/const DEBUG_VERBOSE = false/' "$DEP_FILE"
echo "DEBUG_VERBOSE restored to false"
