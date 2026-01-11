#!/bin/bash
# Setup full bioBakery databases for production use
# Location: scripts/setup_full_databases.sh

set -e

DB_DIR="${1:-databases}"

echo "================================================"
echo "Setting up FULL bioBakery databases"
echo "================================================"
echo ""
echo "Target directory: $DB_DIR"
echo ""
echo "WARNING: This will download ~60 GB of data"
echo "Estimated time: 2-4 hours depending on connection"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 0
fi

# Check for container runtime
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "Error: Neither apptainer nor singularity found"
    exit 1
fi

echo "Using container runtime: $CONTAINER_CMD"
echo ""

mkdir -p "$DB_DIR"

# MetaPhlAn database
echo "================================================"
echo "Downloading MetaPhlAn database (~20 GB)..."
echo "================================================"
echo ""
mkdir -p "$DB_DIR/metaphlan_db"
$CONTAINER_CMD exec docker://quay.io/biocontainers/metaphlan:4.2.4--pyhdfd78af_0 \
    metaphlan --install --db_dir "$DB_DIR/metaphlan_db"

# HUMAnN databases
echo ""
echo "================================================"
echo "Downloading HUMAnN ChocoPhlAn database (~17 GB)..."
echo "================================================"
echo ""
mkdir -p "$DB_DIR/humann_db"
$CONTAINER_CMD exec docker://biobakery/humann:latest \
    humann_databases --update-config no --download chocophlan full "$DB_DIR/humann_db"

echo ""
echo "================================================"
echo "Downloading HUMAnN UniRef90 database (~20 GB)..."
echo "================================================"
echo ""
$CONTAINER_CMD exec docker://biobakery/humann:latest \
    humann_databases --update-config no --download uniref uniref90_diamond "$DB_DIR/humann_db"

echo ""
echo "================================================"
echo "Downloading HUMAnN utility mapping (~1.4 GB)..."
echo "================================================"
echo ""
$CONTAINER_CMD exec docker://biobakery/humann:latest \
    humann_databases --update-config no --download utility_mapping full "$DB_DIR/humann_db"

# KneadData database
echo ""
echo "================================================"
echo "Downloading KneadData human genome (~4 GB)..."
echo "================================================"
echo ""
mkdir -p "$DB_DIR/kneaddata_db"
$CONTAINER_CMD exec docker://biobakery/kneaddata:latest \
    kneaddata_database --download human_genome bowtie2 "$DB_DIR/kneaddata_db"

echo ""
echo "================================================"
echo "✓ All databases downloaded successfully!"
echo "================================================"
echo ""
echo "Database sizes:"
du -sh "$DB_DIR"/*
echo ""
echo "Total size:"
du -sh "$DB_DIR"
echo ""
echo "To use in pipeline:"
echo "  nextflow run main.nf \\"
echo "    --kneaddata_db_dir $DB_DIR/kneaddata_db \\"
echo "    --metaphlan_db_dir $DB_DIR/metaphlan_db \\"
echo "    --humann_db_dir $DB_DIR/humann_db"
