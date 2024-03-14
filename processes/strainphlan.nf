#!/usr/bin/env nextflow

process strainphlan_sample2markers {
    tag "strainphlan sample2markers"
    container "$params.metaphlan_image"

    input:
    path metaphlan_db
    path sams

    output:
    path "consensus_markers/*.pkl"

    """
    mkdir -p consensus_markers
    sample2markers.py \
    	-i $sams \
	-d $metaphlan_db \
	-o consensus_markers \
        -n ${task.cpus} \
    """
}

process strainphlan_extractmarkers {
    tag "strainphlan extract_markers"
    container "$params.metaphlan_image"

    input:
    path metaphlan_db
    val clade

    output:
    path("markers/*")

    """
    mkdir -p markers
    extract_markers.py -d ${metaphlan_db}/*.pkl -c ${clade} -o markers/
    """
}

process strainphlan {
    container "${params.metaphlan_image}"
    publishDir "$params.outdir/strainphlan", mode: "copy", overwrite: true

    input:
    path db
    path consensus_markers
    val clade
    path markers
    path reference_genomes

    output:
    path "${clade}/*", emit: all
    tuple val("${clade}"), path("${clade}/*.tre"), emit: tre

    script:
    """
    mkdir -p ${clade}
    if [ -d $reference_genomes ] && (( \$(find $reference_genomes | wc -l) > 0 )); then
        FLAGS="-r ${reference_genomes}/*"
    else
        FLAGS=""
    fi
    
    mkdir -p local_tmp
    
    strainphlan \
        -s ${consensus_markers} \
        -d ${db} \
        -m ${markers} \
        -o ${clade} \
        -n ${task.cpus} \
        -c ${clade} \
        --non_interactive \
        --mutation_rates \
        --tmp local_tmp \
        \$FLAGS \
    """
}
