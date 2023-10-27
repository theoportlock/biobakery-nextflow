#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kneaddata; kneaddata_init; kneaddata_summary } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_init } from './processes/metaphlan.nf'
include { humann; humann_init; humann_regroup; humann_rename } from './processes/humann.nf'

//include { strainphlan; strainphlan_sample2markers; strainphlan_extractmarkers } from './processes/strainphlan.nf'

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.input)

    kneaddata_db = kneaddata_init(params.kneaddata_db)
    metaphlan_db = metaphlan_init(params.metaphlan_db)
    //humann_db = humann_init()
    
    knead_out = kneaddata(read_pairs_ch, kneaddata_db)
    //kneaddata_summary(knead_out.log.collect())
    metaphlan_out = metaphlan(knead_out[0], knead_out[1], metaphlan_db)
    //metaphlan_merge(metaphlan_out.profile.collect())
    //humann_out = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], humann_db)
    //regroup_out = humann_regroup(humann_out[0], humann_out[1])
    //humann_rename(regroup_out)
}
