#!/bin/bash

config='conf/nesi.config'

module load Nextflow/23.10.0

nextflow main.nf \
	-c $config \
	-resume \
	-bg \
	-with-tower \
	--project uoa03941 \
	--input 'testdata/*R{1,2}*' \
	--k \
	--m

