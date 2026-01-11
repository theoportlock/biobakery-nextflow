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
    STRAINPHLAN_EXTRACT_MARKERS;
    STRAINPHLAN_EXTRACT_REF_MARKERS;
    STRAINPHLAN_BUILD_TREES;
    STRAINPHLAN_TRANSMISSION;
    STRAINPHLAN_AGGREGATE_RESULTS
} from './modules/local/strainphlan/main.nf'

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

    // Create variables to store metaphlan outputs and database
    def metaphlan_out
    def metaphlan_db
    def precomputed_profiles

    // Initialize Metaphlan if requested
    if (params.m) {
        metaphlan_db = METAPHLAN_INIT(params.metaphlan_db)
        // If kneaddata was run, merge the cleaned reads; otherwise merge the original pairs.
        merged_reads = params.k ? MERGE(kneaddata_out.cleaned_reads) : MERGE(read_pairs)
        
        metaphlan_out = METAPHLAN_RUN(merged_reads, metaphlan_db)
	
	// Extract profiles and mapouts from metaphlan tuple output
	metaphlan_profiles = metaphlan_out.profile.map { sample, profile, mapout -> profile }
	metaphlan_mapouts = metaphlan_out.profile.map { sample, profile, mapout -> mapout }
	
	METAPHLAN_MERGE(metaphlan_profiles.collect())

        metaphlan_to_gtdb_out = METAPHLAN_TO_GTDB(metaphlan_out.profile)
	metaphlan_profiles_gtdb = metaphlan_to_gtdb_out.profile.map { sample, profile, mapout -> profile }
	METAPHLAN_MERGE_GTDB(metaphlan_profiles_gtdb.collect())
    }

    // Run StrainPhlAn if requested and MetaPhlAn was enabled
    if (params.s && params.m) {
        // Validate required parameters
        if (params.sgb_list == null || params.sgb_list.isEmpty()) {
            error "StrainPhlAn requires --sgb_list to be specified (file with SGB IDs, one per line)"
        }
        if (params.strain_metadata == null || params.strain_metadata.isEmpty()) {
            error "StrainPhlAn requires --strain_metadata to be specified (TSV with sample_id, subject, relation, timepoint columns)"
        }

        // Extract markers from SAM files
        strainphlan_markers = STRAINPHLAN_EXTRACT_MARKERS(metaphlan_out.samout)

        // Extract reference markers if reference genomes are provided
        def ref_markers
        if (params.reference_genomes_dir != null && !params.reference_genomes_dir.isEmpty()) {
            ref_genomes = Channel.fromPath(params.reference_genomes_dir)
            ref_markers = STRAINPHLAN_EXTRACT_REF_MARKERS(ref_genomes).ref_markers
        } else {
            // Create empty channel if no references
            ref_markers = Channel.fromPath('NO_REF')
        }

        // Read SGB list and create channel
        sgb_list = Channel.fromPath(params.sgb_list)
            .splitText()
            .map { it.trim() }
            .filter { it.length() > 0 && !it.startsWith('#') }

        // Collect all markers into a single directory
        all_markers = strainphlan_markers.map { sample, pkl -> pkl }.collect()

        // Build trees for each SGB
        strainphlan_trees = STRAINPHLAN_BUILD_TREES(
            sgb_list,
            all_markers,
            ref_markers.collect()
        )

        // Load metadata
        metadata = Channel.fromPath(params.strain_metadata)

        // Run transmission analysis for each SGB
        strainphlan_transmission = STRAINPHLAN_TRANSMISSION(
            strainphlan_trees.tree,
            metadata.first()
        )

        // Aggregate results
        STRAINPHLAN_AGGREGATE_RESULTS(
            strainphlan_transmission.transmission.map { sgb, info -> info }.collect()
        )
    } else if (params.s && !params.m) {
        error "StrainPhlAn requires MetaPhlAn to be enabled (--m flag must be set)"
    }

    // Load pre-computed metaphlan profiles if humann is requested without metaphlan
    if (params.h && !params.m) {
        if (params.metaphlan_profile_dir == null || params.metaphlan_profile_dir.isEmpty()) {
            error "HUMANN requires either Metaphlan to be enabled (--m) or pre-computed profiles directory (--metaphlan_profile_dir) to be specified."
        }
        if (params.metaphlan_bowtie2_db_dir == null || params.metaphlan_bowtie2_db_dir.isEmpty()) {
            error "HUMANN requires MetaPhlAn bowtie2 database directory (--metaphlan_bowtie2_db_dir) to be specified."
        }

        // Load pre-computed profiles
        precomputed_profiles = Channel.fromPath("${params.metaphlan_profile_dir}/*_kneaddata_paired_profile.tsv")
            .map { profile_file ->
                def sample_id = profile_file.name.replaceAll(/_kneaddata_paired_profile\.tsv$/, '')
                tuple(sample_id, profile_file)
            }

        // Merge reads if kneaddata was run, otherwise use raw reads
        merged_reads = params.k ? MERGE(kneaddata_out.cleaned_reads) : MERGE(read_pairs)
    }

    // Initialize HUMANN if requested
    if (params.h) {
        // Determine metaphlan database and profiles based on whether metaphlan was run
        if (params.m) {
            // Using metaphlan output
            humann_db = HUMANN_INIT()
            
            // Join merged reads with metaphlan profiles by sample ID
            humann_input = merged_reads.join(metaphlan_out.profile)
                .map { sample, reads, profile, mapout -> tuple(sample, reads, profile) }
            
            humann_out = HUMANN_RUN(humann_input, humann_db, metaphlan_db)
        } else {
            // Using pre-computed profiles
            if (params.metaphlan_bowtie2_db_dir == null || params.metaphlan_bowtie2_db_dir.isEmpty()) {
                error "HUMANN requires MetaPhlAn bowtie2 database directory (--metaphlan_bowtie2_db_dir) when using pre-computed profiles."
            }
            
            humann_db = HUMANN_INIT()
            
            // Match precomputed profiles with merged reads by sample ID
            humann_input = merged_reads.join(precomputed_profiles)
                .map { sample, reads, profile -> tuple(sample, reads, profile) }
            
            humann_out = HUMANN_RUN(humann_input, humann_db, params.metaphlan_bowtie2_db_dir, params.metaphlan_bowtie2_db_dir)
        }
        
        humann_merged = HUMANN_MERGE(
            humann_out.genefamilies.map { sample, file -> file }.collect(),
            humann_out.pathabundance.map { sample, file -> file }.collect(),
            humann_out.pathcoverage.map { sample, file -> file }.collect()
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

