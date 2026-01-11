#!/bin/bash

# Minimal test - just MetaPhlAn with test data
# Use this to verify the basic pipeline works before adding StrainPhlAn

set -e

CONFIG='conf/laptop.config'
OUTDIR='test_output_basic'

echo "================================================"
echo "Running basic MetaPhlAn test (no StrainPhlAn)..."
echo "================================================"

# Clean previous output
[ -d "$OUTDIR" ] && rm -rf "$OUTDIR"

nextflow run main.nf \
    -c $CONFIG \
    -resume \
    --input 'testdata/*_L001_R{1,2}_*.fastq.gz' \
    --m \
    --outdir "$OUTDIR"

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Basic test passed!"
    echo "  Check ${OUTDIR}/metaphlan/ for outputs"
else
    echo "✗ Test failed with exit code: $EXIT_CODE"
fi

exit $EXIT_CODE
