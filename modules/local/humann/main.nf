process HUMANN {
    label "process_medium"
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
    label "error_retry"

    tag "humann - merge outputs"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefamilies
    path pathabundance
    path pathcoverage

    output:
    path "humann_genefamilies_rpk.tsv", emit: genefamilies
    path "humann_pathabundance_rpk.tsv", emit: pathabundance
    path "humann_pathcoverage.tsv", emit: pathcoverage

    script:
    """
    humann_join_tables -i . -o humann_genefamilies_rpk.tsv --file_name genefamilies
    humann_join_tables -i . -o humann_pathabundance_rpk.tsv --file_name pathabundance
    humann_join_tables -i . -o humann_pathcoverage.tsv --file_name pathcoverage
    """
}

process HUMANN_RENORM {
    label "process_medium"
    label "error_retry"

    tag "humann - renormalize to cpm"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefamilies
    path pathabundance

    output:
    path "humann_genefamilies_cpm.tsv", emit: genefamilies
    path "humann_pathabundance_cpm.tsv", emit: pathabundance

    script:
    """
    humann_renorm_table -i $genefamilies -u "cpm" -o humann_genefamilies_cpm.tsv
    humann_renorm_table -i $pathabundance -u "cpm" -o humann_pathabundance_cpm.tsv
    """
}


process HUMANN_REGROUP {
    label "process_medium"
    label "error_retry"

    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefamilies
    path humann_db

    output:
    path "humann_level4ec.tsv", emit: ec
    path "humann_ko.tsv", emit: ko
    path "humann_pfam.tsv", emit: pfam
    path "humann_go.tsv", emit: go
    path "humann_eggnog.tsv", emit: eggnog

    script:
    """
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_level4ec_uniref90.txt.gz --output humann_level4ec.tsv
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_ko_uniref90.txt.gz --output humann_ko.tsv
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_pfam_uniref90.txt.gz --output humann_pfam.tsv
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_go_uniref90.txt.gz --output humann_go.tsv
    humann_regroup_table --input $genefamilies -c ${humann_db}/utility_mapping/map_eggnog_uniref90.txt.gz --output humann_eggnog.tsv
    """
}

process HUMANN_RENAME {
    label "process_medium"
    label "error_retry"

    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann", mode: "copy", overwrite: true
    container "$params.humann_image"

    input:
    path genefamilies
    path pathcoverage
    path pathabundance
    path ec
    path ko
    path pfam
    path go

    output:
    path "humann_genefamilies_rename.tsv"
    path "humann_pathcoverage_rename.tsv"
    path "humann_pathabundance_rename.tsv"
    path "humann_ec_rename.tsv"
    path "humann_ko_rename_orthology.tsv"
    path "humann_ko_rename_pathway.tsv"
    path "humann_ko_rename_module.tsv"
    path "humann_pfam_rename.tsv"
    path "humann_go_rename.tsv"

    script:
    """
    humann_rename_table -i $genefamilies -n metacyc-rxn -s -o humann_genefamilies_rename.tsv
    humann_rename_table -i $pathcoverage -n metacyc-pwy -s -o humann_pathcoverage_rename.tsv
    humann_rename_table -i $pathabundance -n metacyc-pwy -s -o humann_pathabundance_rename.tsv
    humann_rename_table -i $ec -n ec -s -o humann_ec_rename.tsv
    humann_rename_table -i $ko -n kegg-orthology -s -o humann_ko_rename_orthology.tsv
    humann_rename_table -i $ko -n kegg-pathway -s -o humann_ko_rename_pathway.tsv
    humann_rename_table -i $ko -n kegg-module -s -o humann_ko_rename_module.tsv
    humann_rename_table -i $pfam -n pfam -s -o humann_pfam_rename.tsv
    humann_rename_table -i $go -n infogo1000 -s -o humann_go_rename.tsv
    """
}

