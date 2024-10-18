process HUMANN {
    label "process_high"
    label "process_long"
    label "error_retry"

    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    container "$params.humann_image"

    input:
    tuple val(sample), path(reads)
    path profile
    path metaphlan_db
    path humann_db

    output:
    val sample
    path "${sample}_genefamilies.tsv"  , emit: genefamilies
    path "${sample}_pathabundance.tsv" , emit: pathabundance
    path "${sample}_pathcoverage.tsv"  , emit: pathcoverage

    script:
    """
    humann \
        --input ${reads} \
        --taxonomic-profile $profile \
        --threads ${task.cpus} \
        --remove-temp-output \
        --nucleotide-database ${humann_db}/chocophlan \
        --protein-database ${humann_db}/uniref \
	--metaphlan-options "--bowtie2db $metaphlan_db" \
        --output-basename $sample \
        --output .
    """
}
process HUMANN_INIT {
    label "process_single"

    tag "humann download databases"
    container "$params.humann_image"

    output:
    path "humann_db"

    script:
    if (params.humann_db_dir == null || params.humann_db_dir.isEmpty()) {
        """
        humann_databases --update-config no --download chocophlan full humann_db
        humann_databases --update-config no --download uniref uniref90_diamond humann_db
        humann_databases --update-config no --download utility_mapping full humann_db
        """
    } else {
        """
        ln -s ${params.humann_db_dir} humann_db
        """
    }
}

process HUMANN_MERGE {
    label "process_medium"

    tag "humann - merge outputs"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefam
    path pathway
    path pathway_coverage

    output:
    path "humann3_genefam_rpk.tsv", emit: genefam
    path "humann3_pathway_rpk.tsv", emit: pathway
    path "humann3_pathway_coverage.tsv", emit: pathway_coverage

    script:
    """
    humann_join_tables -i genefam -o humann3_genefam_rpk.tsv
    humann_join_tables -i pathway -o  humann3_pathway_rpk.tsv
    humann_join_tables -i pathway_coverage -o  humann3_pathway_coverage.tsv
    """
}

process HUMANN_RENORM {
    label "process_medium"

    tag "humann - renormalize to cpm"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefam
    path pathway

    output:
    path "humann3_genefam_cpm.tsv", emit: genefam
    path "humann3_pathway_cpm.tsv", emit: pathway

    script:
    """
    # Renormalise counts from RPK (reads per kilobase) to CPM (copies per million)
    humann_renorm_table -i humann3_genefam_rpk.tsv -u "cpm" -o humann3_genefam_cpm.tsv
    humann_renorm_table -i humann3_pathway_rpk.tsv -u "cpm" -o humann3_pathway_cpm.tsv
    """
}


process HUMANN_REGROUP {
    label "process_single"

    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"
    container "$params.humann_image"

    input:
    path genefam
    path humann_db

    output:
    path "humann2_ecs.tsv"
    path "humann2_kos.tsv"
    path "humann2_pfam.tsv"

    script:
    """
    humann_regroup_table --input $genefam -c ${humann_db}/utility_mapping/map_level4ec_uniref90.txt.gz --output humann2_ecs.tsv
    humann_regroup_table --input $genefam -c ${humann_db}/utility_mapping/map_ko_uniref90.txt.gz --output humann2_kos.tsv
    humann_regroup_table --input $genefam -c ${humann_db}/utility_mapping/map_pfam_uniref90.txt.gz --output humann2_pfam.tsv

    humann_regroup_table --input demo_fastq/demo_genefamilies-cpm.tsv  --output demo_fastq/rxn-cpm.tsv --groups uniref90_rxn
    """
}

process HUMANN_RENAME {
    label "process_single"

    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"
    container "$params.humann_image"

    input:
    path rxn

    output:
    path "humann3_rxn_rename.tsv"

    script:
    """
    humann_rename_table \
        --names metacyc-rxn \
        --input $rxn \
        --output ${sample}_rxn_rename.tsv
    """
}

