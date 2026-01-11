# Test Data

This directory contains test data for the bioBakery Nextflow pipeline.

## Files

### FASTQ Files
- `FG00004_S26_L001_R{1,2}_001.fastq.gz` - Donor sample, lane 1
- `FG00004_S26_L002_R{1,2}_001.fastq.gz` - Pre-FMT recipient, lane 2  
- `FG00004_S26_L003_R{1,2}_001.fastq.gz` - Post-FMT recipient, lane 3
- `FG00004_S26_L004_R{1,2}_001.fastq.gz` - Additional sample, lane 4
- `FG00005_S50_L001_R{1,2}_001.fastq.gz` - Second donor sample

### StrainPhlAn Configuration

**sgb_list.txt** - Species Genome Bins for strain analysis
```
t__SGB1234
t__SGB5678
t__SGB9012
```

**strain_metadata.tsv** - Sample metadata for transmission analysis
| sample_id | subject | relation | timepoint |
|-----------|---------|----------|-----------|
| FG00004_S26_L001 | donor | donor | 0 |
| FG00004_S26_L002 | recipient | baseline | 0 |
| FG00004_S26_L003 | recipient | followup | 30 |
| FG00005_S50_L001 | donor | donor | 0 |

## Usage

### Basic Test (MetaPhlAn only)
```bash
./test_basic.sh
```

### Full StrainPhlAn Test
```bash
./test_strainphlan.sh
```

### Custom Test
```bash
nextflow run main.nf \
    --input 'testdata/*_L001_R{1,2}_*.fastq.gz' \
    --m --s \
    --sgb_list testdata/sgb_list.txt \
    --strain_metadata testdata/strain_metadata.tsv \
    --outdir test_output
```

## Notes

- The FASTQ files should contain real metagenomic data for meaningful testing
- The SGB IDs in `sgb_list.txt` should match species actually present in your samples
- Sample IDs in `strain_metadata.tsv` must match the FASTQ filename prefixes (e.g., `FG00004_S26_L001` from `FG00004_S26_L001_R1_001.fastq.gz`)
- For StrainPhlAn to work, you need at least 2-3 samples with sufficient sequencing depth
