nextflow.enable.dsl=2

include {
    HUMANN_RUN;
    HUMANN_INIT;
    HUMANN_MERGE;
    HUMANN_RENORM;
    HUMANN_REGROUP;
    HUMANN_RENAME
} from '../../modules/local/humann/main.nf'

workflow test_humann_run {
    take:
    reads_profile  // channel: tuple val(sample), path(reads), path(profile)
    db             // path: humann_db
    metaphlan_db   // val: metaphlan_db name
    bowtie2_db     // val: bowtie2_db_dir
    
    main:
    HUMANN_RUN(reads_profile, db, metaphlan_db, bowtie2_db)
    
    emit:
    genefamilies = HUMANN_RUN.out.genefamilies
    pathabundance = HUMANN_RUN.out.pathabundance
    pathcoverage = HUMANN_RUN.out.pathcoverage
}

workflow test_humann_init {
    main:
    HUMANN_INIT()
    
    emit:
    db = HUMANN_INIT.out
}

workflow test_humann_merge {
    take:
    genefamilies   // channel: path(genefamilies)
    pathabundance  // channel: path(pathabundance)
    pathcoverage   // channel: path(pathcoverage)
    
    main:
    HUMANN_MERGE(genefamilies, pathabundance, pathcoverage)
    
    emit:
    genefamilies = HUMANN_MERGE.out.genefamilies
    pathabundance = HUMANN_MERGE.out.pathabundance
    pathcoverage = HUMANN_MERGE.out.pathcoverage
}

workflow test_humann_renorm {
    take:
    genefamilies   // path: genefamilies
    pathabundance  // path: pathabundance
    
    main:
    HUMANN_RENORM(genefamilies, pathabundance)
    
    emit:
    genefamilies = HUMANN_RENORM.out.genefamilies
    pathabundance = HUMANN_RENORM.out.pathabundance
}

workflow test_humann_regroup {
    take:
    genefamilies  // path: genefamilies
    db            // path: humann_db
    
    main:
    HUMANN_REGROUP(genefamilies, db)
    
    emit:
    ec = HUMANN_REGROUP.out.ec
    ko = HUMANN_REGROUP.out.ko
    pfam = HUMANN_REGROUP.out.pfam
    go = HUMANN_REGROUP.out.go
    eggnog = HUMANN_REGROUP.out.eggnog
}
