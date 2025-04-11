process MERGE {
    label "process_single"

    tag "merge reads for metaphlan and humann"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_merged.fastq.gz")

    script:
    """
    # Concatenate paired-end reads into one file
    cat ${reads[0]} ${reads[1]} > ${sample}_merged.fastq.gz
    """
}
