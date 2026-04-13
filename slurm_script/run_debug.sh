#!/bin/bash

# Run test1-3 with DEBUG_VERBOSE=true for manhattan and euclidean norms, then restore

BASE_DIR="/home/htc/aleoputra/project/FPmeetsFW"
TESTCASE_DIR="$BASE_DIR/testcase"
DEP_FILE="$BASE_DIR/dependencies.jl"
OUT_BASE="$BASE_DIR/testcase/output"

# Clean previous results
rm -rf "$OUT_BASE"
mkdir -p "$OUT_BASE/manhattan" "$OUT_BASE/euclidean"

# Flip DEBUG_VERBOSE to true
sed -i 's/const DEBUG_VERBOSE = false/const DEBUG_VERBOSE = true/' "$DEP_FILE"
echo "DEBUG_VERBOSE set to true"

for NORM in manhattan euclidean; do
    for i in 1 2 3; do
        INSTANCE="$TESTCASE_DIR/test${i}.mps"
        OUTPUT="$OUT_BASE/$NORM/test${i}_output.txt"
        echo "Running test${i} [$NORM]..."
        julia --project "$BASE_DIR/run_test.jl" "$INSTANCE" "$NORM" 0.5 vanilla agnostic > "$OUTPUT" 2>&1
        echo "Done. Output: $OUTPUT"
    done
done

# Restore DEBUG_VERBOSE to false
sed -i 's/const DEBUG_VERBOSE = true/const DEBUG_VERBOSE = false/' "$DEP_FILE"
echo "DEBUG_VERBOSE restored to false"
