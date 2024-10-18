process KNEADDATA {
    label "process_medium"
    label "error_retry"

    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"
    container "$params.kneaddata_image"

    input:
    tuple val(sample), path(reads)
    path kneaddata_db

    output:
    tuple val(sample), path("${sample}_kneaddata_paired_{1,2}.fastq.gz") , emit: cleaned_reads
    path "${sample}_kneaddata.log"                                       , emit: log

    shell:
    """
    kneaddata \
    	--input ${reads[0]} \
	--input ${reads[1]} \
	--reference-db $kneaddata_db \
	--processes ${task.cpus} \
	--output-prefix ${sample}_kneaddata \
	--output .
    gzip *.fastq
    """  
}

process KNEADDATA_INIT {
    label "process_single"

    tag "kneaddata build database"
    container "$params.kneaddata_image"

    input:
    val kneaddata_db

    output:
    path "kneaddata_db"

    script:
    if (params.kneaddata_db_dir == null || params.kneaddata_db_dir.isEmpty()) {
        """
        echo "Kneaddata database directory not found, downloading database..."
        kneaddata_database --download $kneaddata_db bowtie2 kneaddata_db
        """
    } else {
        """
        echo "Using existing Kneaddata database directory: ${params.kneaddata_db_dir}"
        ln -s ${params.kneaddata_db_dir} kneaddata_db
        """
    }
}

process KNEADDATA_SUMMARY {
    label "process_single"

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
    kneaddata_read_count_table
    	--input knead_log_dir
	--output kneaddata_summary.tsv
    """  
}
