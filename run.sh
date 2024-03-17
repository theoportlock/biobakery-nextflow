#!/bin/bash
module load Nextflow Singularity

PROJECT=''
FILEREGEX=''
nextflow ~/biobakery-nextflow/main.nf \
        --project $PROJECT \
	--input $FILEREGEX \
        -c nesi.config \
        -with-tower \
	-resume \
        -bg
