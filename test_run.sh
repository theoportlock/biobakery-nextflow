#!/bin/bash
nextflow main.nf \
	-c conf/laptop.config \
	-resume \
	-bg \
	--input 'biobakery-nextflow/testdata/*R{1,2}*' \
	--k \
	--m \
	--h

