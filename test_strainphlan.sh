#!/bin/bash

# Test script for StrainPhlAn integration
# This script runs the pipeline with the test data in testdata/

set -e

CONFIG='conf/laptop.config'
OUTDIR='test_output_strainphlan'

echo "================================================"
echo "Running StrainPhlAn integration test..."
echo "================================================"
echo ""
echo "Test data:"
echo "  - FASTQ files: testdata/FG00004_S26_L00*_R{1,2}_001.fastq.gz"
echo "  - SGB list: testdata/sgb_list.txt"
echo "  - Metadata: testdata/strain_metadata.tsv"
echo ""
echo "Output directory: ${OUTDIR}"
echo ""

# Clean previous test output
if [ -d "$OUTDIR" ]; then
    echo "Cleaning previous test output..."
    rm -rf "$OUTDIR"
fi

echo "Starting pipeline..."
echo ""

nextflow run main.nf \
    -c $CONFIG \
    -resume \
    --input 'testdata/*_L001_R{1,2}_*.fastq.gz' \
    --m \
    --s \
    --sgb_list testdata/sgb_list.txt \
    --strain_metadata testdata/strain_metadata.tsv \
    --strain_threshold 0.1 \
    --min_samples_per_marker 1 \
    --min_markers_per_sample 5 \
    --save_strain_distances \
    --outdir "$OUTDIR"

EXIT_CODE=$?

echo ""
echo "================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Test completed successfully!"
    echo ""
    echo "Expected outputs:"
    echo "  - ${OUTDIR}/metaphlan/*.sam.bz2 (SAM files for marker extraction)"
    echo "  - ${OUTDIR}/strainphlan/markers/*.pkl (extracted markers)"
    echo "  - ${OUTDIR}/strainphlan/t__SGB*/ (per-SGB results)"
    echo "  - ${OUTDIR}/strainphlan/strainphlan_transmission_summary.csv"
else
    echo "✗ Test failed with exit code: $EXIT_CODE"
    echo ""
    echo "Check .nextflow.log for details"
fi
echo "================================================"

exit $EXIT_CODE
