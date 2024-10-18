#!/bin/bash
module load Nextflow Apptainer

nextflow main.nf \
	-c conf/nesi.config \
	-with-tower \
	-resume \
	-w /nesi/nobackup/uoa03902/biobakery4/work \
	-bg \
	--input '/scale_wlg_persistent/filesets/home/tpor598/biobakery-nextflow/testdata/*R{1,2}*' \
	--outdir /nesi/nobackup/uoa03902/biobakery4/output \
        --project uoa03624 \
	--k \
	--m \
	--h

