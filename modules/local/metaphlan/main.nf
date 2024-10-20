process METAPHLAN {
    label "process_medium"

    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv}"
    container "$params.metaphlan_image"

    input:
    tuple val(sample), path(reads)
    path metaphlan_db

    output:
    val sample, emit: sample
    path "${sample}_profile.tsv", emit: profile
    path "${sample}.sam.bz2"

    script:
    """
    metaphlan \
        --bowtie2out ${sample}_bowtie2.tsv \
        --samout ${sample}.sam.bz2 \
        --input_type fastq \
        --nproc ${task.cpus} \
        --bowtie2db $metaphlan_db \
        --index ${params.metaphlan_db} \
        ${reads} \
        ${sample}_profile.tsv
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
    if (params.metaphlan_db_dir == null || params.metaphlan_db_dir.isEmpty()) {
        """
	metaphlan \
	    --install \
	    --index ${metaphlan_db} \
	    --bowtie2db metaphlan_db
        """
    } else {
        """
        ln -s ${params.metaphlan_db_dir} metaphlan_db
        """
    }
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

process METAPHLAN_MERGE_GTDB {
    label "process_single"

    tag "metaphlan merge outputs"
    container "$params.metaphlan_image"
    publishDir "$params.outdir/metaphlan", mode: "copy", overwrite: true

    input:
    path metaphlan_profiles

    output:
    path "metaphlan_merged_profiles_gtdb.tsv"

    script:
    """
    merge_metaphlan_tables.py $metaphlan_profiles > metaphlan_merged_profiles_gtdb.tsv
    """
}

process METAPHLAN_TO_GTDB {
    label "process_single"

    tag "convert metaphlan sgb to gtdb"
    container "$params.metaphlan_image"

    input:
    val(sample)
    path(profile)

    output:
    path "${sample}_profile_gtdb.tsv", emit: profile

    script:
    """
    sgb_to_gtdb_profile.py \
        -i $profile \
	-d ${params.metaphlan_db}.tsv \
	-o ${sample}_profile_gtdb.tsv
    """
}
