#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { KNEADDATA; KNEADDATA_INIT; KNEADDATA_SUMMARY } from './modules/local/kneaddata/main.nf'
include { METAPHLAN; METAPHLAN_INIT; METAPHLAN_MERGE} from './modules/local/metaphlan/main.nf'
include { HUMANN; HUMANN_INIT; HUMANN_REGROUP; HUMANN_RENAME; HUMANN_MERGE } from './modules/local/humann/main.nf'

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.input)

    if( params.kneaddata_db_dir == null ) {
        kneaddata_db = KNEADDATA_INIT(params.kneaddata_db)
        } else {
        kneaddata_db = file(params.kneaddata_db_dir)
        }
    knead_out = KNEADDATA(read_pairs_ch, kneaddata_db)
    KNEADDATA_SUMMARY(knead_out.log.collect())

    //if( params.metaphlan_db_dir == null ) {
    //    metaphlan_db = METAPHLAN_INIT(params.metaphlan_db)
    //    } else {
    //    metaphlan_db = file(params.metaphlan_db_dir)
    //    }
    //metaphlan_out = METAPHLAN(knead_out[0], knead_out[1], metaphlan_db)
    //METAPHLAN_MERGE(metaphlan_out.profile.collect())

    ////humann_db = HUMANN_INIT(params.humann_db_dir)
    //humann_db = HUMANN_INIT()
    //humann_out = HUMANN(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], metaphlan_db, humann_db)
    ////regroup_out = HUMANN_REGROUP(metaphlan_out[0], humann_out[0], humann_out[1])
    ////humann_rename(regroup_out)
    //HUMANN_MERGE(humann_out[1].collect(),humann_out[2].collect(),humann_out[3].collect())
}
