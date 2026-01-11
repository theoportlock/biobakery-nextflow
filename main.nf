#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { 
    KNEADDATA_RUN; 
    KNEADDATA_INIT; 
    KNEADDATA_SUMMARY 
} from './modules/local/kneaddata/main.nf'

include { 
    METAPHLAN_RUN; 
    METAPHLAN_INIT; 
    METAPHLAN_MERGE; 
    METAPHLAN_TO_GTDB; 
    METAPHLAN_MERGE_GTDB 
} from './modules/local/metaphlan/main.nf'

include { 
    HUMANN_RUN; 
    HUMANN_INIT; 
    HUMANN_MERGE; 
    HUMANN_RENORM; 
    HUMANN_REGROUP; 
    HUMANN_RENAME 
} from './modules/local/humann/main.nf'

include { 
    MERGE 
} from './modules/local/mergereads/main.nf'

workflow {

    // Read input file pairs from the parameter
    read_pairs = Channel.fromFilePairs(params.input)

    // Initialize kneaddata if requested
    if (params.k) {
        kneaddata_db  = KNEADDATA_INIT(params.kneaddata_db)
        kneaddata_out = KNEADDATA_RUN(read_pairs, kneaddata_db)
        KNEADDATA_SUMMARY(kneaddata_out.log.collect())
    }

    // Create a variable for merged_reads that will be used by both Metaphlan and Humann
    def merged_reads

    // Initialize Metaphlan if requested
    if (params.m) {
        metaphlan_db = METAPHLAN_INIT(params.metaphlan_db)
        // If kneaddata was run, merge the cleaned reads; otherwise merge the original pairs.
        merged_reads = params.k ? MERGE(kneaddata_out.cleaned_reads) : MERGE(read_pairs)
        
        metaphlan_out = METAPHLAN_RUN(merged_reads, metaphlan_db)
	metaphlan_profiles_list = metaphlan_out.profile
		.map { it.toAbsolutePath().toString() }
		.collectFile(name: 'metaphlan_abs_paths.list', newLine: true, sort: true)
	METAPHLAN_MERGE(metaphlan_profiles_list)

        metaphlan_to_gtdb_out = METAPHLAN_TO_GTDB(metaphlan_out.sample, metaphlan_out.profile)
	metaphlan_profiles_gtdb_list = metaphlan_to_gtdb_out.profile
		.map { it.toAbsolutePath().toString() }
		.collectFile(name: 'metaphlan_gtdb_abs_paths.list', newLine: true, sort: true)
	METAPHLAN_MERGE_GTDB(metaphlan_profiles_gtdb_list)
    }

    // Initialize HUMANN if requested – note HUMANN requires Metaphlan output
    if (params.h) {
        // Enforce that HUMANN cannot run without Metaphlan
        if ( !params.m ) {
            error "HUMANN requires Metaphlan to be enabled. Please use '--m' along with '--h'."
        }

        humann_db = HUMANN_INIT()
        humann_out = HUMANN_RUN(merged_reads, metaphlan_out.profile, metaphlan_db, humann_db)
        
        humann_merged = HUMANN_MERGE(
            humann_out.genefamilies.collect(),
            humann_out.pathabundance.collect(),
            humann_out.pathcoverage.collect()
        )
        humann_renorm  = HUMANN_RENORM(humann_merged.genefamilies, humann_merged.pathabundance)
        humann_regroup = HUMANN_REGROUP(humann_renorm.genefamilies, humann_db)
        
        HUMANN_RENAME(
            humann_renorm.genefamilies,
            humann_merged.pathcoverage,
            humann_renorm.pathabundance,
            humann_regroup.ec,
            humann_regroup.ko,
            humann_regroup.pfam,
            humann_regroup.go
        )
    }
}

