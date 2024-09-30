#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { KNEADDATA; KNEADDATA_INIT; KNEADDATA_SUMMARY } from './modules/local/kneaddata/main.nf'
include { METAPHLAN; METAPHLAN_INIT; METAPHLAN_MERGE } from './modules/local/metaphlan/main.nf'
include { HUMANN; HUMANN_INIT; HUMANN_MERGE } from './modules/local/humann/main.nf'

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.input)

    // Initialize kneaddata
    if (params.k && params.kneaddata_db_dir == null) {
        kneaddata_db = KNEADDATA_INIT(params.kneaddata_db)
    } else if (params.k) {
        kneaddata_db = file(params.kneaddata_db_dir)
    }

    if (params.k) {
        knead_out = KNEADDATA(read_pairs_ch, kneaddata_db)
        KNEADDATA_SUMMARY(knead_out.log.collect())
    }

    // Initialize metaphlan
    if (params.m && params.metaphlan_db_dir == null) {
        metaphlan_db = METAPHLAN_INIT(params.metaphlan_db)
    } else if (params.m) {
        metaphlan_db = file(params.metaphlan_db_dir)
    }

    if (params.m && params.k) {
        metaphlan_out = METAPHLAN(knead_out[0], knead_out[1], metaphlan_db)
        METAPHLAN_MERGE(metaphlan_out.profile.collect())
    } else if (params.m) {
        metaphlan_out = METAPHLAN(read_pairs_ch, read_pairs_ch, metaphlan_db) // Run directly on read_pairs if kneaddata isn't used
        METAPHLAN_MERGE(metaphlan_out.profile.collect())
    }

    // Initialize humann
    if (params.h && params.m) {
        humann_db = HUMANN_INIT()
        humann_out = HUMANN(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], metaphlan_db, humann_db)
        HUMANN_MERGE(humann_out[1].collect(), humann_out[2].collect(), humann_out[3].collect())
    } else if (params.h) {
        humann_db = HUMANN_INIT()
        humann_out = HUMANN(read_pairs_ch, read_pairs_ch, read_pairs_ch, file(''), humann_db) // Run directly on input if kneaddata/metaphlan aren't used
        HUMANN_MERGE(humann_out[1].collect(), humann_out[2].collect(), humann_out[3].collect())
    }
}
