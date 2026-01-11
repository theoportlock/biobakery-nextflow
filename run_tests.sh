#!/bin/bash
# BioBakery Pipeline Test Runner
# Runs all nf-test suites with prerequisite checks
# Location: run_tests.sh

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Default options
RUN_MODULE_TESTS=true
RUN_WORKFLOW_TESTS=true
RUN_INTEGRATION_TESTS=true
VERBOSE=false
CLEAN_AFTER=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --module-only)
            RUN_WORKFLOW_TESTS=false
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        --workflow-only)
            RUN_MODULE_TESTS=false
            RUN_INTEGRATION_TESTS=false
            shift
            ;;
        --integration-only)
            RUN_MODULE_TESTS=false
            RUN_WORKFLOW_TESTS=false
            RUN_INTEGRATION_TESTS=true
            shift
            ;;
        --verbose|-v)
            VERBOSE=true
            shift
            ;;
        --clean)
            CLEAN_AFTER=true
            shift
            ;;
        --help|-h)
            echo "BioBakery Pipeline Test Runner"
            echo ""
            echo "Usage: ./run_tests.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --module-only         Run only module unit tests"
            echo "  --workflow-only       Run only standalone workflow tests"
            echo "  --integration-only    Run only integration workflow tests"
            echo "  --verbose, -v         Show detailed test output"
            echo "  --clean               Remove test outputs after completion"
            echo "  --help, -h            Show this help message"
            echo ""
            echo "Examples:"
            echo "  ./run_tests.sh                    # Run all tests"
            echo "  ./run_tests.sh --module-only      # Run only module tests"
            echo "  ./run_tests.sh --verbose --clean  # Run all tests with cleanup"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run './run_tests.sh --help' for usage information"
            exit 1
            ;;
    esac
done

echo "================================================"
echo "BioBakery Pipeline Test Runner"
echo "================================================"
echo ""

# Check if nf-test is installed
if ! command -v nf-test &> /dev/null; then
    echo -e "${RED}✗ nf-test not found${NC}"
    echo ""
    echo "Please install nf-test:"
    echo "  curl -fsSL https://code.askimed.com/install/nf-test | bash"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓${NC} nf-test is installed"

# Check database prerequisites
echo ""
echo "Checking database prerequisites..."
echo ""

DB_ISSUES=0

# Check KneadData database
if [ ! -d "testdata/kneaddata_db" ] || [ ! -f "testdata/kneaddata_db/test_host.1.bt2" ]; then
    echo -e "${YELLOW}⚠${NC} KneadData database not found"
    echo "  Run: ./scripts/setup_kneaddata_db.sh"
    DB_ISSUES=$((DB_ISSUES + 1))
else
    echo -e "${GREEN}✓${NC} KneadData database found"
fi

# Check MetaPhlAn database
if [ ! -d "testdata/metaphlan_db" ]; then
    echo -e "${YELLOW}⚠${NC} MetaPhlAn database not found"
    echo "  This database is large (~20 GB)"
    echo "  Either run: ./scripts/setup_full_databases.sh"
    echo "  Or use existing database with --metaphlan_db_dir parameter"
    DB_ISSUES=$((DB_ISSUES + 1))
else
    echo -e "${GREEN}✓${NC} MetaPhlAn database found"
fi

# Check HUMAnN demo database
if [ ! -d "testdata/humann_demo_db" ] || [ ! -d "testdata/humann_demo_db/chocophlan" ]; then
    echo -e "${YELLOW}⚠${NC} HUMAnN demo database not found"
    echo "  Run: ./scripts/setup_demo_databases.sh"
    DB_ISSUES=$((DB_ISSUES + 1))
else
    echo -e "${GREEN}✓${NC} HUMAnN demo database found"
fi

# Check precomputed MetaPhlAn outputs
if [ ! -d "testdata/precomputed_metaphlan/metaphlan" ]; then
    echo -e "${YELLOW}⚠${NC} Pre-computed MetaPhlAn outputs not found"
    echo "  Run: ./scripts/setup_precomputed_metaphlan.sh"
    echo "  (Optional: only needed for HUMAnN standalone tests)"
    # Don't increment DB_ISSUES for this one as it's optional
else
    echo -e "${GREEN}✓${NC} Pre-computed MetaPhlAn outputs found"
fi

if [ $DB_ISSUES -gt 0 ]; then
    echo ""
    echo -e "${RED}✗ $DB_ISSUES required database(s) missing${NC}"
    echo ""
    echo "Some tests will fail without these databases."
    echo "Set up databases with scripts in scripts/ directory."
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo ""
echo "================================================"
echo "Running Tests"
echo "================================================"
echo ""

VERBOSE_FLAG=""
if [ "$VERBOSE" = true ]; then
    VERBOSE_FLAG="--verbose"
fi

TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to run tests and track results
run_test() {
    local test_file=$1
    local test_name=$2
    
    echo "Running: $test_name"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    if nf-test test "$test_file" $VERBOSE_FLAG; then
        echo -e "${GREEN}✓${NC} PASSED: $test_name"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo -e "${RED}✗${NC} FAILED: $test_name"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    echo ""
}

# Run module unit tests
if [ "$RUN_MODULE_TESTS" = true ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Module Unit Tests"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    run_test "tests/modules/mergereads.nf.test" "MergeReads Module"
    run_test "tests/modules/kneaddata.nf.test" "KneadData Module"
    run_test "tests/modules/metaphlan.nf.test" "MetaPhlAn Module"
    run_test "tests/modules/humann.nf.test" "HUMAnN Module"
    run_test "tests/modules/strainphlan.nf.test" "StrainPhlAn Module"
fi

# Run standalone workflow tests
if [ "$RUN_WORKFLOW_TESTS" = true ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Standalone Workflow Tests"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    run_test "tests/workflow/kneaddata_workflow.nf.test" "KneadData Workflow (--k)"
    run_test "tests/workflow/metaphlan_workflow.nf.test" "MetaPhlAn Workflow (--m)"
    run_test "tests/workflow/humann_workflow.nf.test" "HUMAnN Workflow (--h)"
fi

# Run integration workflow tests
if [ "$RUN_INTEGRATION_TESTS" = true ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Integration Workflow Tests"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    run_test "tests/workflow/kneaddata_metaphlan.nf.test" "K+M Integration"
    run_test "tests/workflow/metaphlan_humann.nf.test" "M+H Integration"
    run_test "tests/workflow/strainphlan_workflow.nf.test" "M+S Integration"
    run_test "tests/workflow/full_pipeline_kmh.nf.test" "K+M+H Full Pipeline"
    run_test "tests/workflow/full_pipeline_kms.nf.test" "K+M+S Full Pipeline"
    run_test "tests/workflow/full_pipeline_all.nf.test" "K+M+H+S Complete Pipeline"
fi

# Clean up if requested
if [ "$CLEAN_AFTER" = true ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Cleaning up test outputs..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    
    rm -rf .nf-test/
    rm -rf work/
    rm -rf .nextflow*
    
    echo -e "${GREEN}✓${NC} Test outputs cleaned"
    echo ""
fi

# Print summary
echo "================================================"
echo "Test Summary"
echo "================================================"
echo ""
echo "Total tests run:  $TOTAL_TESTS"
echo -e "Passed:          ${GREEN}$PASSED_TESTS${NC}"
echo -e "Failed:          ${RED}$FAILED_TESTS${NC}"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}✗ $FAILED_TESTS test(s) failed${NC}"
    echo ""
    echo "Review the output above for details on failed tests."
    exit 1
fi
