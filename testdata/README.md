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
## Database Directories (Not in Git)

The following database directories must be created using setup scripts:

### `kneaddata_db/`
Minimal test database for host contamination removal (~10 KB)

**Setup:**
```bash
./scripts/setup_kneaddata_db.sh
```

**Contents:**
- `test_host.fa` - Synthetic 10kb genome
- `test_host.*.bt2` - Bowtie2 index files

### `metaphlan_db/`
MetaPhlAn marker database (~20 GB)

**Setup:**
```bash
./scripts/setup_full_databases.sh
```

Or use existing database via `--metaphlan_db_dir` parameter.

**Contents:**
- `mpa_vJan25_CHOCOPhlAnSGB_202503.*.bt2` - Bowtie2 index files
- `mpa_vJan25_CHOCOPhlAnSGB_202503.pkl` - Marker metadata
- Other supporting files

### `humann_demo_db/`
HUMAnN DEMO databases for testing (~1.5 GB)

**Setup:**
```bash
./scripts/setup_demo_databases.sh
```

**Contents:**
- `chocophlan/` - Demo pangenome database
- `uniref/` - Demo protein database (Diamond format)
- `utility_mapping/` - Gene to pathway mappings (required)

### `precomputed_metaphlan/`
Pre-computed MetaPhlAn outputs for faster testing (optional)

**Setup:**
```bash
./scripts/setup_precomputed_metaphlan.sh /path/to/metaphlan_db
```

**Contents:**
- `metaphlan/` - MetaPhlAn output directory
    - `*_profile.tsv` - Taxonomic profiles
    - `*.sam.bz2` - SAM alignment files (for StrainPhlAn)
    - `*.mapout.bz2` - Bowtie2 mapping info

**Purpose:** Allows testing HUMAnN and StrainPhlAn without re-running MetaPhlAn each time.

## Database Sizes

| Database | Size | Required For | Setup Time |
|----------|------|--------------|------------|
| kneaddata_db | ~10 KB | KneadData tests | < 1 min |
| metaphlan_db | ~20 GB | MetaPhlAn, StrainPhlAn tests | 1-2 hours |
| humann_demo_db | ~1.5 GB | HUMAnN tests | 5-10 min |
| precomputed_metaphlan | ~500 MB | HUMAnN standalone tests (optional) | 10-30 min |

## Using Custom Databases

Instead of creating databases in `testdata/`, you can point to existing databases:

```bash
# In test files or command line
nextflow run main.nf \
    --kneaddata_db_dir /path/to/kneaddata/db \
    --metaphlan_db_dir /path/to/metaphlan/db \
    --humann_db_dir /path/to/humann/db \
    ...
```

## Notes

- Database directories are **not tracked in git** (see `.gitignore`)
- Large files: Consider using symlinks to shared database locations
- For CI/CD: Database artifacts should be cached between runs
- DEMO databases (HUMAnN) are sufficient for testing but not for production
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
