#!/bin/bash
#module load Nextflow Singularity

nextflow main.nf \
	-c laptop.config \
	--input testdata/rawfastq/*R{1,2}*.fastq.gz \
	--k
