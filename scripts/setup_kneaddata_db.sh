#!/bin/bash
# Setup minimal KneadData database for testing
# Creates a small synthetic genome with bowtie2 index
# Location: scripts/setup_kneaddata_db.sh

set -e

# Configuration
DB_DIR="testdata/kneaddata_db"
GENOME_NAME="test_host"
GENOME_SIZE=10000  # 10kb synthetic genome

echo "================================================"
echo "Setting up KneadData test database"
echo "================================================"
echo ""
echo "Target directory: ${DB_DIR}"
echo "Genome size: ${GENOME_SIZE} bp"
echo "Estimated time: 1-2 minutes"
echo ""

# Create directory
mkdir -p "$DB_DIR"

# Check for required tools
if ! command -v bowtie2-build &> /dev/null; then
    echo "Error: bowtie2-build not found"
    echo "Please install bowtie2 to continue"
    echo ""
    echo "Install with:"
    echo "  conda install -c bioconda bowtie2"
    echo "  or"
    echo "  sudo apt-get install bowtie2"
    exit 1
fi

echo "Creating synthetic host genome (${GENOME_SIZE} bp)..."

# Generate synthetic genome with Python
python3 - <<EOF
import random
random.seed(42)  # Reproducible genome

with open('${DB_DIR}/${GENOME_NAME}.fa', 'w') as f:
    f.write(">test_host_chr1\n")
    seq = ''.join(random.choices('ACGT', k=${GENOME_SIZE}))
    # Write in 80-char lines (FASTA format)
    for i in range(0, len(seq), 80):
        f.write(seq[i:i+80] + '\n')

print(f"✓ Synthetic genome created: ${DB_DIR}/${GENOME_NAME}.fa")
EOF

echo ""
echo "Building bowtie2 index..."
bowtie2-build -q "${DB_DIR}/${GENOME_NAME}.fa" "${DB_DIR}/${GENOME_NAME}"

echo ""
echo "================================================"
echo "Validating database installation..."
echo "================================================"
echo ""

# Validate bowtie2 index
VALID=true

if [ ! -f "${DB_DIR}/${GENOME_NAME}.fa" ]; then
    echo "✗ FASTA file not found"
    VALID=false
else
    echo "✓ FASTA file exists"
fi

# Check for all bowtie2 index files
for ext in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
    if [ ! -f "${DB_DIR}/${GENOME_NAME}.${ext}" ]; then
        echo "✗ Bowtie2 index file ${GENOME_NAME}.${ext} not found"
        VALID=false
    fi
done

if [ "$VALID" = true ]; then
    echo "✓ All bowtie2 index files present"
    
    echo ""
    echo "================================================"
    echo "✓ KneadData test database created successfully!"
    echo "================================================"
    echo ""
    echo "Database files:"
    ls -lh "$DB_DIR"
    echo ""
    echo "To use in pipeline:"
    echo "  nextflow run main.nf --k --kneaddata_db_dir $DB_DIR --kneaddata_db $GENOME_NAME"
    echo ""
    echo "To use in tests:"
    echo "  nf-test test tests/workflow/kneaddata_workflow.nf.test"
else
    echo ""
    echo "✗ Database validation failed"
    echo "Please check the errors above and try again"
    exit 1
fi
