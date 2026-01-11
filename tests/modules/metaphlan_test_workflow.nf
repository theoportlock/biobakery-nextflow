nextflow.enable.dsl=2

include { 
    METAPHLAN_RUN; 
    METAPHLAN_INIT; 
    METAPHLAN_MERGE;
    METAPHLAN_TO_GTDB;
    METAPHLAN_MERGE_GTDB
} from '../../modules/local/metaphlan/main.nf'

workflow test_metaphlan_run {
    take:
    reads  // channel: tuple val(sample), path(reads)
    db     // path: metaphlan_db
    
    main:
    METAPHLAN_RUN(reads, db)
    
    emit:
    profile = METAPHLAN_RUN.out.profile
    samout = METAPHLAN_RUN.out.samout
}

workflow test_metaphlan_init {
    take:
    db_name  // val: database name
    
    main:
    METAPHLAN_INIT(db_name)
    
    emit:
    db = METAPHLAN_INIT.out
}

workflow test_metaphlan_merge {
    take:
    profiles  // channel: path(profiles)
    
    main:
    METAPHLAN_MERGE(profiles)
    
    emit:
    merged = METAPHLAN_MERGE.out
}

workflow test_metaphlan_to_gtdb {
    take:
    profile  // channel: tuple val(sample), path(profile), path(mapout)
    
    main:
    METAPHLAN_TO_GTDB(profile)
    
    emit:
    profile = METAPHLAN_TO_GTDB.out.profile
}

workflow test_metaphlan_merge_gtdb {
    take:
    profiles  // channel: path(profiles)
    
    main:
    METAPHLAN_MERGE_GTDB(profiles)
    
    emit:
    merged = METAPHLAN_MERGE_GTDB.out
}
