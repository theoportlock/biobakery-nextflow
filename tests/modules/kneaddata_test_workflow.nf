nextflow.enable.dsl=2

include { 
    KNEADDATA_RUN; 
    KNEADDATA_INIT; 
    KNEADDATA_SUMMARY 
} from '../../modules/local/kneaddata/main.nf'

workflow test_kneaddata_run {
    take:
    reads  // channel: tuple val(sample), path(reads)
    db     // path: kneaddata_db
    
    main:
    KNEADDATA_RUN(reads, db)
    
    emit:
    cleaned_reads = KNEADDATA_RUN.out.cleaned_reads
    log = KNEADDATA_RUN.out.log
}

workflow test_kneaddata_init {
    take:
    db_name  // val: database name
    
    main:
    KNEADDATA_INIT(db_name)
    
    emit:
    db = KNEADDATA_INIT.out
}

workflow test_kneaddata_summary {
    take:
    logs  // channel: path(logs)
    
    main:
    KNEADDATA_SUMMARY(logs)
    
    emit:
    summary = KNEADDATA_SUMMARY.out
}
