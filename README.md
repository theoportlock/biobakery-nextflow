# BioBakery Nextflow Pipeline
Author: Theo Portlock, based on projects from Kevin Bonham and Tommi Vatanen

This Nextflow pipeline is designed for metagenomic data processing using tools from the [bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/Home) suite: **KneadData**, **MetaPhlAn**, **HUMAnN**, and **StrainPhlAn**. It enables modular execution of these tools with initialization, processing, and post-processing steps for microbiome profiling and strain-level tracking.

## Features

- Conditional execution using flags: `--k` for KneadData, `--m` for MetaPhlAn, `--h` for HUMAnN, and `--s` for StrainPhlAn.
- Automatic database initialization.
- Output merging, renormalization, regrouping, and renaming for downstream analysis.
- Integration with GTDB.
- Strain-level tracking and transmission analysis for FMT studies.

## Usage

### Basic Usage (All Tools)
```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --k --m --h
```

### KneadData + MetaPhlAn Only
```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --k --m
```

### StrainPhlAn for Strain Tracking
StrainPhlAn requires MetaPhlAn to be enabled (`--m`) and additional input files:

```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --m --s \
    --sgb_list sgb_list.txt \
    --strain_metadata strain_metadata.tsv \
    --strain_threshold 0.1
```

**Required parameters:**
- `--sgb_list`: Text file with SGB IDs (one per line, e.g., `t__SGB1234`)
- `--strain_metadata`: TSV file with columns: `sample_id`, `subject`, `relation`, `timepoint`

**Optional parameters:**
- `--reference_genomes_dir`: Directory with reference genome `.fna` files
- `--strain_threshold`: Distance threshold for transmission events (default: 0.1)
- `--min_samples_per_marker`: Minimum samples a marker must appear in (default: 2)
- `--min_markers_per_sample`: Minimum markers required per sample (default: 10)
- `--save_strain_distances`: Save pairwise distance matrices

**Example files:**

`sgb_list.txt`:
```
t__SGB1234
t__SGB5678
t__SGB9012
```

`strain_metadata.tsv`:
```
sample_id	subject	relation	timepoint
donor_1	donor	donor	0
pre_fmt_1	recipient	baseline	0
post_fmt_1	recipient	followup	30
```

**Outputs:**
- Per-SGB phylogenetic trees (`RAxML_bestTree.*.tre`)
- Multiple sequence alignments (`.aln` files)
- Transmission events per SGB (`transmission_events.info`)
- Aggregated summary across all SGBs (`strainphlan_transmission_summary.csv`)

### HUMAnN with Pre-computed MetaPhlAn Profiles
If you have already run MetaPhlAn separately and have pre-computed taxonomic profiles, you can run HUMAnN without re-running MetaPhlAn:

```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --h \
    --metaphlan_profile_dir ./metaphlan \
    --metaphlan_bowtie2_db_dir /path/to/metaphlan/bowtie2db
```

**Expected profile file naming:** `{sample_id}_kneaddata_paired_profile.tsv`  
For example: `B1-01_kneaddata_paired_profile.tsv`, `B1-02_kneaddata_paired_profile.tsv`, etc.

**Database parameters:**
- `--metaphlan_profile_dir`: Directory containing pre-computed MetaPhlAn profiles
- `--metaphlan_bowtie2_db_dir`: Path to MetaPhlAn bowtie2 database (e.g., `mpa_vJan25_CHOCOPhlAnSGB_202503`)

### Custom Database Directories
For any tool, you can specify custom database paths:
```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --k --m --h \
    --kneaddata_db_dir /custom/kneaddata/db \
    --metaphlan_db_dir /custom/metaphlan/db \
    --humann_db_dir /custom/humann/db
```

## Pipeline Modules

### KneadData (`--k`)
- `KNEADDATA_INIT`: Initializes KneadData with the provided database.
- `KNEADDATA`: Runs KneadData on paired reads.
- `KNEADDATA_SUMMARY`: Summarizes KneadData logs.

### MetaPhlAn (`--m`)
- `METAPHLAN_INIT`: Initializes MetaPhlAn with the provided database.
- `MERGE`: Merges cleaned (or raw) reads.
- `METAPHLAN_RUN`: Runs MetaPhlAn on merged reads, generating profiles and SAM files.
- `METAPHLAN_MERGE`: Merges MetaPhlAn profiles.
- `METAPHLAN_TO_GTDB`: Maps MetaPhlAn profiles to GTDB taxonomy.
- `METAPHLAN_MERGE_GTDB`: Merges GTDB-mapped MetaPhlAn profiles.

### StrainPhlAn (`--s`)
- Requires `--m` to be enabled.
- `STRAINPHLAN_EXTRACT_MARKERS`: Extracts strain markers from SAM files.
- `STRAINPHLAN_EXTRACT_REF_MARKERS`: Extracts markers from reference genomes (optional).
- `STRAINPHLAN_BUILD_TREES`: Builds phylogenetic trees for each SGB.
- `STRAINPHLAN_TRANSMISSION`: Analyzes transmission events using distance thresholds.
- `STRAINPHLAN_AGGREGATE_RESULTS`: Aggregates transmission events across all SGBs.

### HUMAnN (`--h`)
- Can run with MetaPhlAn enabled (`--m`) or with pre-computed profiles (`--metaphlan_profile_dir`).
- `HUMANN_INIT`: Initializes HUMAnN with required databases.
- `HUMANN_RUN`: Runs HUMAnN using MetaPhlAn output (when `--m` is enabled).
- `HUMANN_RUN_PRECOMPUTED`: Runs HUMAnN using pre-computed MetaPhlAn profiles (when `--metaphlan_profile_dir` is provided).
- `HUMANN_MERGE`: Merges HUMAnN output files.
- `HUMANN_RENORM`: Renormalizes merged output.
- `HUMANN_REGROUP`: Regroups gene families into various annotation schemes (EC, KO, PFAM, GO).
- `HUMANN_RENAME`: Renames and organizes final outputs.

## Notes

- Input must be in paired format, e.g., `sample_1.fastq.gz` and `sample_2.fastq.gz`.

## License

MIT License
