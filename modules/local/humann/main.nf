process HUMANN {
    label "process_medium"

    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    container "$params.humann_image"

    input:
    val sample
    path profile
    path catkneads
    path metaphlan_db
    path humann_db

    output:
    val sample
    path "${sample}_genefamilies.tsv"
    path "${sample}_pathabundance.tsv"
    path "${sample}_pathcoverage.tsv"

    script:
    """
    humann --input $catkneads \
        --taxonomic-profile $profile \
        --output . \
        --threads ${task.cpus} \
        --remove-temp-output \
        --nucleotide-database ${humann_db}/chocophlan \
        --protein-database ${humann_db}/uniref \
	--metaphlan-options "--bowtie2db $metaphlan_db" \
        --output-basename $sample
    """
}

process HUMANN_INIT {
    label "process_single"

    tag "humann download databases"
    container "$params.humann_image"

    output:
    path "humann_db"

    script:
    """
    humann_databases --update-config no --download chocophlan full humann_db
    humann_databases --update-config no --download uniref uniref90_diamond humann_db
    humann_databases --update-config no --download utility_mapping full humann_db
    """
}   

process HUMANN_REGROUP {
    label "process_single"

    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"
    container "$params.humann_image"

    input:
    val sample
    path genefamilies
    path humann_db

    output:
    val sample, emit: sample
    path "${sample}_ecs.tsv"
    path "${sample}_kos.tsv"
    path "${sample}_pfams.tsv"

    script:
    """
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_level4ec_uniref90.txt.gz --output ${sample}_ecs.tsv 
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_ko_uniref90.txt.gz --output ${sample}_kos.tsv 
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_pfam_uniref90.txt.gz --output ${sample}_pfams.tsv 
    """
}   

process HUMANN_RENAME {
    label "process_single"

    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"
    container "$params.humann_image"

    input:
    val sample
    path rxn

    output:
    val sample, emit: sample
    path "${sample}_rxn_rename.tsv"

    script:
    """
    humann_rename_table --input $rxn --output ${sample}_rxn_rename.tsv --names metacyc-rxn
    """
}

process HUMANN_MERGE {
    label "process_medium"

    tag "humann_merge"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefamilies
    path pathabundance
    path pathcoverage

    output:
    path "merged_genefamilies.tsv"
    path "merged_pathabundance.tsv"
    path "merged_pathcoverage.tsv"

    script:
    """
    mkdir genefamilies pathabundance pathcoverage
    cp -a $genefamilies genefamilies
    cp -a $pathabundance pathabundance
    cp -a $pathcoverage pathcoverage
    humann_join_tables --input pathabundance --output merged_pathabundance.tsv
    humann_join_tables --input pathcoverage --output merged_pathcoverage.tsv
    humann_join_tables --input genefamilies --output merged_genefamilies.tsv
    """
}
