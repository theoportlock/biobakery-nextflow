# BioBakery Nextflow Pipeline
Author: Theo Portlock, based on projects from Kevin Bonham and Tommi Vatanen

This Nextflow pipeline is designed for metagenomic data processing using tools from the [bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/Home) suite: **KneadData**, **MetaPhlAn**, and **HUMAnN**. It enables modular execution of these tools with initialization, processing, and post-processing steps for microbiome profiling.

## Features

- Conditional execution using flags: `--k` for KneadData, `--m` for MetaPhlAn, and `--h` for HUMAnN.
- Automatic database initialization.
- Output merging, renormalization, regrouping, and renaming for downstream analysis.
- Integration with GTDB.

## Usage

```bash
nextflow run main.nf --input 'data/*_{1,2}.fastq.gz' --k --m --h --kneaddata_db /path/to/kneaddata_db --metaphlan_db /path/to/metaphlan_db
```

## Pipeline Modules

### KneadData (`--k`)
- `KNEADDATA_INIT`: Initializes KneadData with the provided database.
- `KNEADDATA`: Runs KneadData on paired reads.
- `KNEADDATA_SUMMARY`: Summarizes KneadData logs.

### MetaPhlAn (`--m`)
- `METAPHLAN_INIT`: Initializes MetaPhlAn with the provided database.
- `MERGE`: Merges cleaned (or raw) reads.
- `METAPHLAN`: Runs MetaPhlAn on merged reads.
- `METAPHLAN_MERGE`: Merges MetaPhlAn profiles.
- `METAPHLAN_TO_GTDB`: Maps MetaPhlAn profiles to GTDB taxonomy.
- `METAPHLAN_MERGE_GTDB`: Merges GTDB-mapped MetaPhlAn profiles.

### HUMAnN (`--h`)
- Requires MetaPhlAn (`--m`) to be enabled.
- `HUMANN_INIT`: Initializes HUMAnN with required databases.
- `HUMANN`: Runs HUMAnN on merged reads and MetaPhlAn output.
- `HUMANN_MERGE`: Merges HUMAnN output files.
- `HUMANN_RENORM`: Renormalizes merged output.
- `HUMANN_REGROUP`: Regroups gene families into various annotation schemes (EC, KO, PFAM, GO).
- `HUMANN_RENAME`: Renames and organizes final outputs.

## Notes

- Input must be in paired format, e.g., `sample_1.fastq.gz` and `sample_2.fastq.gz`.

## License

MIT License
