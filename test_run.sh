#!/bin/bash
nextflow main.nf \
	-c conf/laptop.config \
	-with-tower \
	-resume \
	--input 'testdata/*R{1,2}*' \
	--k \
	--m

