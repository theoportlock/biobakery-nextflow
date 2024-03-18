process METAPHLAN {
    label "process_medium" 

    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv}"
    container "$params.metaphlan_image"

    input:
    tuple val(sample), path(kneads)
    path unmatched
    path metaphlan_db

    output:
    val sample, emit: sample
    path "${sample}_profile.tsv", emit: profile
    path "${sample}_grouped.fastq.gz"
    path "${sample}_bowtie2.tsv"
    path "${sample}.sam.bz2"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]
    def unf = unmatched[0]
    def unr = unmatched[1]

    """
    cat $forward $reverse $unf $unr > ${sample}_grouped.fastq.gz
    metaphlan ${sample}_grouped.fastq.gz ${sample}_profile.tsv \
        --bowtie2out ${sample}_bowtie2.tsv \
        --samout ${sample}.sam.bz2 \
        --input_type fastq \
        --nproc ${task.cpus} \
        --bowtie2db $metaphlan_db \
	--index ${params.metaphlan_db}
    """
}
 
process METAPHLAN_INIT {
    label "process_single"

    tag "metaphlan install database"
    container "$params.metaphlan_image"

    input:
    val metaphlan_db

    output:
    path "metaphlan_db"

    script:
    """
    metaphlan --install --index ${metaphlan_db} --bowtie2db metaphlan_db
    """
}

process METAPHLAN_MERGE {
    label "process_single"

    tag "metaphlan merge outputs"
    container "$params.metaphlan_image"
    publishDir "$params.outdir/metaphlan", mode: "copy", overwrite: true

    input:
    path metaphlan_profiles

    output:
    path "metaphlan_merged_profiles.tsv"

    script:
    """
    merge_metaphlan_tables.py $metaphlan_profiles > metaphlan_merged_profiles.tsv
    """
}
