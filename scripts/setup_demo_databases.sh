#!/bin/bash
# Setup HUMAnN demo databases for testing
# Location: scripts/setup_demo_databases.sh

set -e

# Configuration
DB_DIR="testdata/humann_demo_db"
CONTAINER_IMAGE="docker://biobakery/humann:latest"

echo "================================================"
echo "Setting up HUMAnN demo databases"
echo "================================================"
echo ""
echo "Target directory: ${DB_DIR}"
echo "Estimated size: ~1.5 GB"
echo "Estimated time: 5-10 minutes"
echo ""

# Create directory
mkdir -p "$DB_DIR"

# Check if apptainer/singularity is available
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "Error: Neither apptainer nor singularity found"
    echo "Please install apptainer or singularity to continue"
    exit 1
fi

echo "Using container runtime: $CONTAINER_CMD"
echo ""

# Download databases using container
echo "Downloading ChocoPhlAn DEMO database..."
$CONTAINER_CMD exec "$CONTAINER_IMAGE" \
    humann_databases --update-config no --download chocophlan DEMO "$DB_DIR"

echo ""
echo "Downloading UniRef90 DEMO database..."
$CONTAINER_CMD exec "$CONTAINER_IMAGE" \
    humann_databases --update-config no --download uniref DEMO_diamond "$DB_DIR"

echo ""
echo "Downloading utility mapping databases..."
$CONTAINER_CMD exec "$CONTAINER_IMAGE" \
    humann_databases --update-config no --download utility_mapping full "$DB_DIR"

echo ""
echo "Note: HUMAnN's utility_mapping has no DEMO variant;"
echo "      the required mapping is the 'full' set (~0.9 GB)."
echo "      This is expected even when using demo ChocoPhlAn/UniRef."

echo ""
echo "================================================"
echo "Validating database installation..."
echo "================================================"
echo ""

# Validate databases
VALID=true

if [ ! -d "$DB_DIR/chocophlan" ]; then
    echo "✗ ChocoPhlAn directory not found"
    VALID=false
else
    echo "✓ ChocoPhlAn directory exists"
fi

if [ ! -d "$DB_DIR/uniref" ]; then
    echo "✗ UniRef directory not found"
    VALID=false
else
    echo "✓ UniRef directory exists"
fi

if [ ! -d "$DB_DIR/utility_mapping" ]; then
    echo "✗ Utility mapping directory not found"
    VALID=false
else
    echo "✓ Utility mapping directory exists"
fi

# Check for specific files
if [ ! -f "$DB_DIR/uniref"/*.dmnd ]; then
    echo "✗ UniRef Diamond database file not found"
    VALID=false
else
    echo "✓ UniRef Diamond database file exists"
fi

if [ "$VALID" = true ]; then
    echo ""
    echo "================================================"
    echo "✓ Demo databases downloaded successfully!"
    echo "================================================"
    echo ""
    echo "Database structure:"
    du -sh "$DB_DIR"/*
    echo ""
    echo "To use in pipeline:"
    echo "  nextflow run main.nf --h --humann_db_dir $DB_DIR"
    echo ""
    echo "To use in tests:"
    echo "  nf-test test tests/workflow/humann_workflow.nf.test"
else
    echo ""
    echo "✗ Database validation failed"
    echo "Please check the errors above and try again"
    exit 1
fi
