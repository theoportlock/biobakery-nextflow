#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kneaddata; kneaddata_init; kneaddata_summary } from './processes/kneaddata.nf'
include { metaphlan; metaphlan_init; metaphlan_merge} from './processes/metaphlan.nf'
include { humann; humann_init; humann_regroup; humann_rename; humann_merge } from './processes/humann.nf'

include { strainphlan_sample2markers } from './processes/strainphlan.nf'
include { strainphlan_extractmarkers as strainphlan_extractmarkersb } from './processes/strainphlan.nf'
include { strainphlan_extractmarkers as strainphlan_extractmarkersl } from './processes/strainphlan.nf'
include { strainphlan as strainphlanb } from './processes/strainphlan.nf'
include { strainphlan as strainphlanl } from './processes/strainphlan.nf'

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.input)

    kneaddata_db = kneaddata_init(params.kneaddata_db)
    metaphlan_db = metaphlan_init(params.metaphlan_db)
    humann_db = humann_init()
    
    knead_out = kneaddata(read_pairs_ch, kneaddata_db)
    kneaddata_summary(knead_out.log.collect())

    metaphlan_out = metaphlan(knead_out[0], knead_out[1], metaphlan_db)
    metaphlan_merge(metaphlan_out.profile.collect())

    strainphlan_markers = strainphlan_sample2markers(metaphlan_db, metaphlan_out[4]).collect()
    strainphlan_markers_output = strainphlan_markers.collect()

    extractedmarkersb = strainphlan_extractmarkersb(metaphlan_db, 't__SGB17278')
    extractedmarkersl = strainphlan_extractmarkersl(metaphlan_db, 't__SGB7144')

    strainphlanb(metaphlan_db, strainphlan_markers, 't__SGB17278', extractedmarkersb, '/scale_wlg_nobackup/filesets/nobackup/uoa03941/Public_genomes/banimalis')
    strainphlanl(metaphlan_db, strainphlan_markers, 't__SGB7144', extractedmarkersl, '/scale_wlg_nobackup/filesets/nobackup/uoa03941/Public_genomes/lrhamnosus')

    humann_out = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], metaphlan_db, humann_db)
    //regroup_out = humann_regroup(metaphlan_out[0], humann_out[0], humann_out[1])
    //humann_rename(regroup_out)
    humann_merge(humann_out[1].collect(),humann_out[2].collect(),humann_out[3].collect())
}
