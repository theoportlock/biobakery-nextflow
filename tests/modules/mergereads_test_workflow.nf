nextflow.enable.dsl=2

include { MERGE } from '../../modules/local/mergereads/main.nf'

workflow test_merge {
    take:
    reads  // channel: tuple val(sample), path(reads)
    
    main:
    MERGE(reads)
    
    emit:
    merged = MERGE.out
}
