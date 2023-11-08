process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    container "$params.humann_image"

    input:
    val sample
    path profile
    path catkneads
    path metaphlan_db
    path humann_nucleotide_db
    path humann_protein_db

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

process humann_init {
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

process humann_regroup {
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

process humann_rename {
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"

    container "biobakery/humann:3.6"

    input:
    val sample
    path rxn

    output:
    val  sample , emit: sample
    path "${sample}_rxn_rename.tsv"

    script:
    """
    humann_rename_table --input $rxn --output ${sample}_rxn_rename.tsv --names metacyc-rxn
    """
}

process humann_merge {
    tag "humann_merge"
    publishDir "$params.outdir/humann"

    container "biobakery/humann:3.6"

    input:
    path humann_profiles

    output:
    val  sample , emit: sample
    path "${sample}_rxn_rename.tsv"

    script:
    """
    group='pathabundances'
    mkdir $group
    ln -s *"$group"* $group
    humann_join_tables --input $group --output humann_merged_${group}.tsv
    """
}
