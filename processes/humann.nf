process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    container "$params.humann_image"

    input:
    val sample
    path profile
    path catkneads
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
        --nucleotide-database $humann_nucleotide_db \
        --protein-database $humann_protein_db \
        --output-basename $sample
    """
}

process humann_init {
    tag "humann download databases"
    container "$params.humann_image"

    output:
    path "humann_nucleotide_db"
    path "humann_protein_db"

    script:
    """
    humann_databases --update-config no --download chocophlan full humann_nucleotide_db
    humann_databases --update-config no --download uniref uniref90_diamond humann_protein_db
    """
}   

process humann_regroup {
    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"
    container "$params.humann_image"

    input:
    val sample
    path genefamilies

    output:
    val sample, emit: sample
    path "${sample}_rxn.tsv"

    script:
    """
    humann_regroup_table --input $genefamilies --output ${sample}_rxn.tsv --groups uniref90_rxn
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
