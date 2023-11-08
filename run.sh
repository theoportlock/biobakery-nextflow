#!/bin/bash
export NXF_SINGULARITY_CACHEDIR=$INSERTDIR
export SINGULARITY_CACHEDIR=$INSERTDIR

module load Nextflow Singularity
nextflow ~/biobakery-nextflow/main.nf \
        --project $INSERTPROJECT \
	-profile engaging \
        -c nesi.config \
        -c run.config \
        -with-tower \
	-resume \
        -bg
