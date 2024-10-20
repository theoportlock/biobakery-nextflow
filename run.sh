#!/bin/bash
module load Nextflow Apptainer

dbdir='/scale_wlg_nobackup/filesets/nobackup/uoa03902/biobakery_db'
nextflow main.nf \
	-c conf/nesi.config \
	-with-tower \
	-resume \
	-bg \
	--input '/scale_wlg_persistent/filesets/home/tpor598/biobakery-nextflow/testdata/*R{1,2}*' \
        --project uoa03624 \
	--kneaddata_db_dir ${dbdir}/kneaddata_db \
	--metaphlan_db_dir ${dbdir}/metaphlan_db \
	--humann_db_dir ${dbdir}/humann_db \
	--k \
	--m \
	--h

