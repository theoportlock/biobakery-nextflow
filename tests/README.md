# StrainPhlAn Testing

This directory contains tests for the StrainPhlAn module and workflow integration.

## Test Structure

- `modules/strainphlan.nf.test` - Unit tests for individual StrainPhlAn processes
- `workflow/strainphlan_workflow.nf.test` - Integration tests for the full workflow

## Running Tests

### Using nf-test framework

```bash
# Install nf-test (if not already installed)
curl -fsSL https://code.askimed.com/install/nf-test | bash

# Run all StrainPhlAn tests
nf-test test tests/modules/strainphlan.nf.test
nf-test test tests/workflow/strainphlan_workflow.nf.test

# Run specific test
nf-test test tests/modules/strainphlan.nf.test --testcase "STRAINPHLAN_EXTRACT_MARKERS"
```

### Using the test script

```bash
# Simple integration test
./test_strainphlan.sh
```

## Test Data Requirements

For the tests to run, you need minimal test data in `testdata/`:

1. **Sample FASTQ files**: `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`, etc.
2. **SGB list**: `sgb_list.txt` (already provided)
3. **Metadata**: `strain_metadata.tsv` (already provided)
4. **Optional reference genomes**: `reference_genomes/*.fna`

## Test Cases

### Module Tests
- Marker extraction from SAM files
- Phylogenetic tree building
- Transmission analysis
- Results aggregation

### Workflow Tests
- Full pipeline: MetaPhlAn → StrainPhlAn
- With reference genomes
- Error handling (missing parameters)
- Error handling (wrong dependencies)

## Expected Outputs

After successful test run:
```
test_output_strainphlan/
├── metaphlan/
│   ├── sample1_profile.tsv
│   └── sample1.sam.bz2
├── strainphlan/
│   ├── markers/
│   │   └── sample1.pkl
│   ├── t__SGB1234/
│   │   ├── RAxML_bestTree.t__SGB1234.StrainPhlAn4.tre
│   │   └── transmission_events.info
│   └── strainphlan_transmission_summary.csv
```

## Cleaning Up

```bash
# Remove test outputs
rm -rf test_output_strainphlan/
rm -rf .nextflow*
rm -rf work/
```
