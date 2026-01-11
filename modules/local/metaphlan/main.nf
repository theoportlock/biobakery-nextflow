process METAPHLAN_RUN {
    label "process_medium"

    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv,*.mapout.bz2,*.sam.bz2}", mode: "copy"
    container "$params.metaphlan_image"

    input:
    tuple val(sample), path(reads)
    path metaphlan_db

    output:
    tuple val(sample), path("${sample}_profile.tsv"), path("${sample}.mapout.bz2"), emit: profile
    tuple val(sample), path("${sample}.sam.bz2"), emit: samout

    script:
    """
    metaphlan \
        ${reads} \
        --input_type fastq \
        --db_dir ${metaphlan_db} \
        --index ${params.metaphlan_db} \
        --bowtie2out ${sample}.mapout.bz2 \
        --samout ${sample}.sam.bz2 \
        --nproc ${task.cpus} \
        -o ${sample}_profile.tsv
    """
}

process METAPHLAN_INIT {
    label "process_single"
    label "process_extra_long"

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
	    --db_dir metaphlan_db
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
    path profiles

    output:
    path "metaphlan_merged_profiles.tsv"

    script:
    """
    merge_metaphlan_tables.py \
        -i ${profiles.join(',')} \
        -o metaphlan_merged_profiles.tsv
    """
}

process METAPHLAN_MERGE_GTDB {
    label "process_single"
    tag "metaphlan merge (GTDB)"

    container "$params.metaphlan_image"
    publishDir "$params.outdir/metaphlan", mode: "copy", overwrite: true

    input:
    path profiles

    output:
    path "metaphlan_merged_profiles_gtdb.tsv"

    script:
    """
    merge_metaphlan_tables.py \
        --gtdb_profiles \
        -i ${profiles.join(',')} \
        -o metaphlan_merged_profiles_gtdb.tsv
    """
}

process METAPHLAN_TO_GTDB {
    label "process_single"

    tag "convert metaphlan sgb to gtdb"
    container "$params.metaphlan_image"

    input:
    tuple val(sample), path(profile), path(mapout)

    output:
    tuple val(sample), path("${sample}_profile_gtdb.tsv"), path(mapout), emit: profile

    script:
    """
    sgb_to_gtdb_profile.py \
        -i $profile \
	-o ${sample}_profile_gtdb.tsv
    """
}
