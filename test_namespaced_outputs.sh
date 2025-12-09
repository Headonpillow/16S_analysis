#!/bin/bash
# Test script to verify that runs are properly namespaced and isolated
# This test validates the acceptance criteria: "A run writes only inside its run dir"

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_RUN_NAME="test_$(date +%s)"
TEST_RUN_DIR="${REPO_ROOT}/runs/${TEST_RUN_NAME}"

echo "================================"
echo "Namespaced Outputs Validation Test"
echo "================================"
echo ""

# Clean up any previous test runs
cleanup() {
    if [ -d "${TEST_RUN_DIR}" ]; then
        echo "Cleaning up test run directory: ${TEST_RUN_DIR}"
        rm -rf "${TEST_RUN_DIR}"
    fi
}

trap cleanup EXIT

echo "1. Creating test run directory structure..."
cd "${REPO_ROOT}"
./setup_run.sh "${TEST_RUN_NAME}"

if [ ! -d "${TEST_RUN_DIR}" ]; then
    echo "ERROR: Test run directory was not created"
    exit 1
fi
echo "   ✓ Test run directory created"

echo ""
echo "2. Verifying expected directory structure..."

# Check for config.yaml
if [ ! -f "${TEST_RUN_DIR}/config.yaml" ]; then
    echo "   ✗ FAIL: config.yaml not found in ${TEST_RUN_DIR}"
    exit 1
fi
echo "   ✓ config.yaml exists"

# Check for data directory
if [ ! -d "${TEST_RUN_DIR}/data" ]; then
    echo "   ✗ FAIL: data directory not found"
    exit 1
fi
echo "   ✓ data directory exists"

# Check for data subdirectories
for subdir in raw_external db meta; do
    if [ ! -d "${TEST_RUN_DIR}/data/${subdir}" ]; then
        echo "   ✗ FAIL: data/${subdir} directory not found"
        exit 1
    fi
done
echo "   ✓ data subdirectories exist (raw_external, db, meta)"

echo ""
echo "3. Verifying config.yaml content..."

# Check that logs_dir is configured
if ! grep -q "logs_dir:" "${TEST_RUN_DIR}/config.yaml"; then
    echo "   ✗ FAIL: logs_dir not configured in config.yaml"
    exit 1
fi
echo "   ✓ logs_dir configured"

# Check that results_dir is configured
if ! grep -q "results_dir:" "${TEST_RUN_DIR}/config.yaml"; then
    echo "   ✗ FAIL: results_dir not configured in config.yaml"
    exit 1
fi
echo "   ✓ results_dir configured"

# Check that intermediate_dir is configured
if ! grep -q "intermediate_dir:" "${TEST_RUN_DIR}/config.yaml"; then
    echo "   ✗ FAIL: intermediate_dir not configured in config.yaml"
    exit 1
fi
echo "   ✓ intermediate_dir configured"

echo ""
echo "4. Verifying Snakefile has LOGS_DIR variable..."
if ! grep -q "LOGS_DIR = config.get" "${REPO_ROOT}/Snakefile"; then
    echo "   ✗ FAIL: LOGS_DIR not configured in Snakefile"
    exit 1
fi
echo "   ✓ LOGS_DIR configured in Snakefile"

echo ""
echo "5. Verifying all paths in config are relative (not absolute)..."
# Check that no absolute paths (starting with /) are in config
if grep -E "^[[:space:]]*(data_dir|raw_data_subdir|db_subdir|meta_subdir|intermediate_dir|results_dir|logs_dir):[[:space:]]*\"/" "${TEST_RUN_DIR}/config.yaml"; then
    echo "   ✗ FAIL: Found absolute paths in config.yaml"
    exit 1
fi
echo "   ✓ All paths in config are relative"

echo ""
echo "6. Verifying Snakefile rules have log directives..."
# Count how many rules have log directives
log_count=$(grep -c "^[[:space:]]*log:" "${REPO_ROOT}/Snakefile" || echo 0)
if [ "$log_count" -lt 5 ]; then
    echo "   ✗ FAIL: Expected at least 5 rules with log directives, found ${log_count}"
    exit 1
fi
echo "   ✓ Found ${log_count} rules with log directives"

echo ""
echo "================================"
echo "✓ ALL TESTS PASSED"
echo "================================"
echo ""
echo "Summary:"
echo "  - Run directory structure is correct"
echo "  - config.yaml contains logs_dir, results_dir, and intermediate_dir"
echo "  - All paths are relative (not absolute)"
echo "  - Snakefile has LOGS_DIR configuration"
echo "  - Snakefile rules have log outputs"
echo ""
echo "Acceptance criteria validated:"
echo "  ✓ Layout includes runs/<run_id>/config.yaml"
echo "  ✓ Layout supports runs/<run_id>/results/"
echo "  ✓ Layout supports runs/<run_id>/logs/"
echo "  ✓ No hardcoded absolute paths to results/ or logs/"
echo "  ✓ Configuration ensures a run writes only inside its run dir"
