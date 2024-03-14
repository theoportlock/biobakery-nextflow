process kneaddata {
    tag "kneaddata $sample"
    //publishDir "$params.outdir/kneaddata"
    //container "https://depot.galaxyproject.org/singularity/kneaddata%3A0.12.0--pyhdfd78af_1"
    container "$params.kneaddata_image"

    input:
    tuple val(sample), path(reads)
    path kneaddata_db

    output:
    tuple val(sample), path("${sample}_kneaddata_paired_{1,2}.fastq.gz")
    path "${sample}_kneaddata_unmatched_{1,2}.fastq.gz"
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    shell:
    """
    kneaddata --input ${reads[0]} --input ${reads[1]} --reference-db $kneaddata_db --output . --processes ${task.cpus} --output-prefix ${sample}_kneaddata
    gzip *.fastq
    """  
}

process kneaddata_init {
    tag "kneaddata build database"
    container "$params.kneaddata_image"

    input:
    val kneaddata_db

    output:
    path "kneaddata_db"

    script:
    """
    kneaddata_database --download $kneaddata_db bowtie2 kneaddata_db/
    """  
}

process kneaddata_summary {
    tag "kneaddata summarise run"
    container "$params.kneaddata_image"
    publishDir "$params.outdir/kneaddata", mode: "copy", overwrite: true

    input:
    path kneaddata_logs

    output:
    path "kneaddata_summary.tsv"

    script:
    """
    mkdir knead_log_dir
    cp $kneaddata_logs knead_log_dir
    kneaddata_read_count_table --input knead_log_dir --output kneaddata_summary.tsv
    """  
}
