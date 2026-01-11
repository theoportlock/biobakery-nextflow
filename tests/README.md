# BioBakery Pipeline Testing

This directory contains comprehensive tests for all BioBakery pipeline modules and workflow integrations.

## Test Structure

### Module Unit Tests (`modules/`)

Unit tests for individual processes in each module:

- `mergereads.nf.test` - Tests for read merging (MERGE process)
- `kneaddata.nf.test` - Tests for host contamination removal (KNEADDATA_RUN, KNEADDATA_SUMMARY, KNEADDATA_INIT)
- `metaphlan.nf.test` - Tests for taxonomic profiling (METAPHLAN_RUN, METAPHLAN_MERGE, METAPHLAN_TO_GTDB, METAPHLAN_INIT)
- `humann.nf.test` - Tests for functional profiling (HUMANN_RUN, HUMANN_MERGE, HUMANN_RENORM, HUMANN_REGROUP)
- `strainphlan.nf.test` - Tests for strain analysis (STRAINPHLAN_EXTRACT_MARKERS, STRAINPHLAN_BUILD_TREES, etc.)

### Workflow Tests (`workflow/`)

#### Standalone Workflow Tests
Tests for individual modules enabled via flags:

- `kneaddata_workflow.nf.test` - KneadData only (`--k`)
- `metaphlan_workflow.nf.test` - MetaPhlAn only (`--m`)
- `humann_workflow.nf.test` - HUMAnN only (`--h` with precomputed profiles)
- `strainphlan_workflow.nf.test` - StrainPhlAn only (`--m --s`)

#### Integration Workflow Tests
Tests for combinations of modules:

- `kneaddata_metaphlan.nf.test` - K+M: Read cleaning → taxonomic profiling
- `metaphlan_humann.nf.test` - M+H: Taxonomic → functional profiling
- `full_pipeline_kmh.nf.test` - K+M+H: Complete functional pipeline
- `full_pipeline_kms.nf.test` - K+M+S: Complete strain pipeline
- `full_pipeline_all.nf.test` - K+M+H+S: All modules enabled

## Prerequisites

Before running tests, you must set up minimal test databases:

### Required Databases

1. **KneadData test database** (~10 KB)
   ```bash
   ./scripts/setup_kneaddata_db.sh
   ```
   Creates a minimal synthetic genome in `testdata/kneaddata_db/`

2. **MetaPhlAn database** (~20 GB)
   ```bash
   ./scripts/setup_full_databases.sh
   ```
   Or use existing database by setting `--metaphlan_db_dir` in tests

3. **HUMAnN DEMO databases** (~1.5 GB)
   ```bash
   ./scripts/setup_demo_databases.sh
   ```
   Downloads minimal ChocoPhlAn, UniRef90, and utility_mapping to `testdata/humann_demo_db/`

4. **Pre-computed MetaPhlAn outputs** (optional, for HUMAnN standalone tests)
   ```bash
   ./scripts/setup_precomputed_metaphlan.sh /path/to/metaphlan_db
   ```
   Generates profiles and SAM files in `testdata/precomputed_metaphlan/`

### Install nf-test

```bash
# Install nf-test (if not already installed)
curl -fsSL https://code.askimed.com/install/nf-test | bash
```

## Running Tests

### Quick Start - Run All Tests

```bash
# Run all tests (module + workflow + integration)
./run_tests.sh
```

### Run Specific Test Categories

```bash
# Run only module unit tests
./run_tests.sh --module-only

# Run only standalone workflow tests
./run_tests.sh --workflow-only

# Run only integration tests
./run_tests.sh --integration-only

# Run with verbose output
./run_tests.sh --verbose

# Run and clean up test outputs after
./run_tests.sh --clean
```

### Run Individual Tests

```bash
# Run specific module test
nf-test test tests/modules/humann.nf.test

# Run specific workflow test
nf-test test tests/workflow/full_pipeline_kmh.nf.test

# Run specific test case
nf-test test tests/modules/metaphlan.nf.test --testcase "METAPHLAN_RUN"
```

### Run Integration Test Scripts (Bash)

```bash
# MetaPhlAn only
./test_basic.sh

# StrainPhlAn integration
./test_strainphlan.sh
```

## Test Data Requirements

### Sample Data (Already Provided)

In `testdata/`:
- **FASTQ files**: Paired-end reads from 5 samples (2 donors, 3 recipient timepoints)
  - `FG00004_S26_L001_R{1,2}_001.fastq.gz` (Donor)
  - `FG00004_S26_L002_R{1,2}_001.fastq.gz` (Pre-FMT recipient)
  - `FG00004_S26_L003_R{1,2}_001.fastq.gz` (Post-FMT recipient)
  - `FG00004_S26_L004_R{1,2}_001.fastq.gz` (Additional sample)
  - `FG00005_S50_L001_R{1,2}_001.fastq.gz` (Second donor)
- **SGB list**: `sgb_list.txt` (Species Genome Bins for StrainPhlAn)
- **Metadata**: `strain_metadata.tsv` (Sample metadata for transmission analysis)

### Databases (Must Be Created)

See **Prerequisites** section above for setup instructions.

## Test Coverage Matrix

| Combination | Test File | Modules Tested | Purpose |
|-------------|-----------|----------------|---------|
| `--k` | kneaddata_workflow.nf.test | KneadData | Quality control |
| `--m` | metaphlan_workflow.nf.test | MetaPhlAn | Taxonomic profiling |
| `--h` | humann_workflow.nf.test | HUMAnN | Functional profiling (with precomputed profiles) |
| `--m --s` | strainphlan_workflow.nf.test | MetaPhlAn + StrainPhlAn | Strain-level analysis |
| `--k --m` | kneaddata_metaphlan.nf.test | KneadData + MetaPhlAn | QC → profiling |
| `--m --h` | metaphlan_humann.nf.test | MetaPhlAn + HUMAnN | Taxonomic → functional |
| `--k --m --h` | full_pipeline_kmh.nf.test | K + M + H | Complete functional pipeline |
| `--k --m --s` | full_pipeline_kms.nf.test | K + M + S | Complete strain pipeline |
| `--k --m --h --s` | full_pipeline_all.nf.test | All modules | Complete bioBakery pipeline |

## Test Categories by Module

### MergeReads Tests
- Read concatenation (R1 + R2 → merged)
- Sample ID preservation
- Output file validation

### KneadData Tests
- Host contamination removal
- Log file generation
- Summary aggregation across samples
- Database initialization

### MetaPhlAn Tests
- Taxonomic profiling from merged reads
- Profile merging across samples
- SGB to GTDB taxonomy conversion
- SAM file generation for StrainPhlAn
- Database initialization

### HUMAnN Tests
- Functional profiling with MetaPhlAn profiles
- Gene family abundance calculation
- Pathway abundance and coverage
- CPM normalization
- Regrouping to functional databases (KO, EC, Pfam, GO, eggNOG)
- Database initialization

### StrainPhlAn Tests
- Marker extraction from SAM files
- Phylogenetic tree construction
- Transmission event analysis
- Results aggregation

## Expected Outputs

### Module Tests
Module tests validate individual process outputs.

### Workflow Tests
After successful workflow test run:
```
test_output/
├── kneaddata/                    # (if --k)
│   ├── *_kneaddata_paired_*.fastq.gz
│   └── kneaddata_summary.tsv
├── metaphlan/
│   ├── sample1_profile.tsv
│   └── sample1.sam.bz2
├── humann/                       # (if --h)
│   ├── humann_genefamilies_cpm.tsv
│   ├── humann_pathabundance_cpm.tsv
│   ├── humann_ko.tsv
│   └── humann_level4ec.tsv
├── strainphlan/
│   ├── markers/
│   │   └── sample1.pkl
│   ├── t__SGB1234/
│   │   ├── RAxML_bestTree.t__SGB1234.StrainPhlAn4.tre
│   │   └── transmission_events.info
│   └── strainphlan_transmission_summary.csv
```

## Troubleshooting

### Test Failures

1. **Database not found errors**
   - Ensure you've run the appropriate setup script from `scripts/`
   - Check that database directories exist in `testdata/`

2. **Container errors**
   - Ensure Apptainer or Singularity is installed
   - Check that container images can be pulled

3. **Memory/resource errors**
   - Tests use config from `conf/laptop.config`
   - Adjust resource limits if needed

4. **Timeout errors**
   - Tests assume databases are pre-setup
   - Database downloads will cause timeouts

### Getting Help

```bash
# Show test runner help
./run_tests.sh --help

# Show nf-test help
nf-test --help

# Run test in verbose mode
nf-test test tests/modules/metaphlan.nf.test --verbose
```

## Cleaning Up

```bash
# Remove test outputs
rm -rf test_output*/
rm -rf .nextflow*
rm -rf work/
rm -rf .nf-test/

# Or use the test runner with --clean flag
./run_tests.sh --clean
```

## CI/CD Integration

The test suite can be integrated into CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
- name: Setup databases
  run: |
    ./scripts/setup_demo_databases.sh
    ./scripts/setup_kneaddata_db.sh

- name: Run tests
  run: ./run_tests.sh --verbose
```

## Contributing

When adding new modules or features:

1. Create module unit tests in `tests/modules/`
2. Add workflow integration tests in `tests/workflow/`
3. Update this README with new test descriptions
4. Ensure all tests pass before submitting PR

Follow the existing test patterns for consistency.
