#!/bin/bash
# Generate pre-computed MetaPhlAn outputs for testing
# This script runs MetaPhlAn once on test samples and saves outputs
# Location: scripts/setup_precomputed_metaphlan.sh

set -e

# Configuration
OUTPUT_DIR="testdata/precomputed_metaphlan"
CONFIG="conf/laptop.config"
METAPHLAN_DB_DIR="${1:-}"  # Optional: path to existing MetaPhlAn database

echo "================================================"
echo "Generating pre-computed MetaPhlAn outputs"
echo "================================================"
echo ""
echo "Target directory: ${OUTPUT_DIR}"
echo "Input data: testdata/*_L001_R{1,2}_*.fastq.gz"
echo ""

if [ -z "$METAPHLAN_DB_DIR" ]; then
    echo "Note: No MetaPhlAn database path provided"
    echo "      Pipeline will download full database (~20 GB, 1-2 hours)"
    echo ""
    read -p "Continue? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 0
    fi
    DB_PARAM=""
else
    echo "Using MetaPhlAn database: $METAPHLAN_DB_DIR"
    DB_PARAM="--metaphlan_db_dir $METAPHLAN_DB_DIR"
fi

echo ""
echo "Running MetaPhlAn on test samples..."
echo "This may take 10-30 minutes depending on database location"
echo ""

# Run MetaPhlAn workflow
nextflow run main.nf \
    -c "$CONFIG" \
    -resume \
    --input 'testdata/*_L001_R{1,2}_*.fastq.gz' \
    --m \
    $DB_PARAM \
    --outdir "$OUTPUT_DIR"

EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo ""
    echo "✗ MetaPhlAn run failed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi

echo ""
echo "================================================"
echo "Validating outputs..."
echo "================================================"
echo ""

# Check for required outputs
VALID=true

if [ ! -d "$OUTPUT_DIR/metaphlan" ]; then
    echo "✗ MetaPhlAn output directory not found"
    VALID=false
else
    echo "✓ MetaPhlAn output directory exists"
fi

# Check for SAM files (required for StrainPhlAn)
SAM_COUNT=$(find "$OUTPUT_DIR/metaphlan" -name "*.sam.bz2" | wc -l)
if [ $SAM_COUNT -eq 0 ]; then
    echo "✗ No SAM files found"
    VALID=false
else
    echo "✓ Found $SAM_COUNT SAM files"
fi

# Check for profile files (required for HUMAnN)
PROFILE_COUNT=$(find "$OUTPUT_DIR/metaphlan" -name "*_profile.tsv" | wc -l)
if [ $PROFILE_COUNT -eq 0 ]; then
    echo "✗ No profile files found"
    VALID=false
else
    echo "✓ Found $PROFILE_COUNT profile files"
fi

if [ "$VALID" = true ]; then
    echo ""
    echo "================================================"
    echo "✓ Pre-computed outputs generated successfully!"
    echo "================================================"
    echo ""
    echo "Output structure:"
    du -sh "$OUTPUT_DIR"/*
    echo ""
    echo "Files generated:"
    find "$OUTPUT_DIR/metaphlan" -type f -name "*.sam.bz2" -o -name "*_profile.tsv" | sort
    echo ""
    echo "Usage:"
    echo "  For HUMAnN tests: --metaphlan_profile_dir $OUTPUT_DIR/metaphlan"
    echo "  For StrainPhlAn: SAM files will be used automatically"
    echo ""
    echo "Note: Add $OUTPUT_DIR to .gitignore (large files)"
else
    echo ""
    echo "✗ Output validation failed"
    echo "Please check the errors above and try again"
    exit 1
fi
