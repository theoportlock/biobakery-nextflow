#!/bin/bash
# Script to update module test infrastructure
set -e

echo "Updating module test infrastructure..."
echo ""

# Create humann_test_workflow.nf
echo "Creating tests/modules/humann_test_workflow.nf..."
cat > tests/modules/humann_test_workflow.nf << 'EOF'
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
EOF

# Update kneaddata.nf.test
echo "Updating tests/modules/kneaddata.nf.test..."
cat > tests/modules/kneaddata.nf.test << 'EOF'
nextflow_workflow {

    name "Test KNEADDATA_RUN Process"
    script "tests/modules/kneaddata_test_workflow.nf"
    workflow "test_kneaddata_run"
    
    test("KNEADDATA_RUN - clean host contamination (Apptainer)") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                kneaddata_image = 'docker://biobakery/kneaddata:latest'
            }
            workflow {
                """
                input[0] = Channel.fromFilePairs('${projectDir}/testdata/*_R{1,2}_001.fastq.gz', size: 2)
                    .first()
                input[1] = file('${projectDir}/testdata/kneaddata_db')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("KNEADDATA_RUN - clean host contamination (Docker)") {
        config "conf/test.docker.config"
        
        when {
            params {
                kneaddata_image = 'docker://biobakery/kneaddata:latest'
            }
            workflow {
                """
                input[0] = Channel.fromFilePairs('${projectDir}/testdata/*_R{1,2}_001.fastq.gz', size: 2)
                    .first()
                input[1] = file('${projectDir}/testdata/kneaddata_db')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test KNEADDATA_INIT Process"
    script "tests/modules/kneaddata_test_workflow.nf"
    workflow "test_kneaddata_init"
    
    test("KNEADDATA_INIT - symlink existing database") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_image = 'docker://biobakery/kneaddata:latest'
            }
            workflow {
                """
                input[0] = 'test_host'
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test KNEADDATA_SUMMARY Process"
    script "tests/modules/kneaddata_test_workflow.nf"
    workflow "test_kneaddata_summary"
    
    test("KNEADDATA_SUMMARY - aggregate log files") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                kneaddata_image = 'docker://biobakery/kneaddata:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/output/kneaddata/*.log')
                    .collect()
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update metaphlan.nf.test
echo "Updating tests/modules/metaphlan.nf.test..."
cat > tests/modules/metaphlan.nf.test << 'EOF'
nextflow_workflow {

    name "Test METAPHLAN_RUN Process"
    script "tests/modules/metaphlan_test_workflow.nf"
    workflow "test_metaphlan_run"
    
    test("METAPHLAN_RUN - profile merged reads (Apptainer)") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
                metaphlan_db = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/*_merged.fastq.gz')
                    .map { file -> tuple(file.simpleName.replaceAll('_merged', ''), file) }
                    .first()
                input[1] = file('${projectDir}/testdata/metaphlan_db')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("METAPHLAN_RUN - profile merged reads (Docker)") {
        config "conf/test.docker.config"
        
        when {
            params {
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
                metaphlan_db = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/*_merged.fastq.gz')
                    .map { file -> tuple(file.simpleName.replaceAll('_merged', ''), file) }
                    .first()
                input[1] = file('${projectDir}/testdata/metaphlan_db')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test METAPHLAN_INIT Process"
    script "tests/modules/metaphlan_test_workflow.nf"
    workflow "test_metaphlan_init"
    
    test("METAPHLAN_INIT - symlink existing database") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
            }
            workflow {
                """
                input[0] = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test METAPHLAN_MERGE Process"
    script "tests/modules/metaphlan_test_workflow.nf"
    workflow "test_metaphlan_merge"
    
    test("METAPHLAN_MERGE - merge multiple profiles") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/precomputed_metaphlan/metaphlan/*_profile.tsv')
                    .collect()
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test METAPHLAN_TO_GTDB Process"
    script "tests/modules/metaphlan_test_workflow.nf"
    workflow "test_metaphlan_to_gtdb"
    
    test("METAPHLAN_TO_GTDB - convert SGB to GTDB") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/precomputed_metaphlan/metaphlan/*_profile.tsv')
                    .map { profile -> 
                        def sample = profile.simpleName.replaceAll('_profile', '')
                        def mapout = file("${projectDir}/testdata/precomputed_metaphlan/metaphlan/\${sample}.mapout.bz2")
                        tuple(sample, profile, mapout)
                    }
                    .first()
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test METAPHLAN_MERGE_GTDB Process"
    script "tests/modules/metaphlan_test_workflow.nf"
    workflow "test_metaphlan_merge_gtdb"
    
    test("METAPHLAN_MERGE_GTDB - merge GTDB profiles") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                metaphlan_image = 'docker://biobakery/metaphlan:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/precomputed_metaphlan/metaphlan/*_profile_gtdb.tsv')
                    .collect()
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update humann.nf.test
echo "Updating tests/modules/humann.nf.test..."
cat > tests/modules/humann.nf.test << 'EOF'
nextflow_workflow {

    name "Test HUMANN_RUN Process"
    script "tests/modules/humann_test_workflow.nf"
    workflow "test_humann_run"
    
    test("HUMANN_RUN - generate functional profiles (Apptainer)") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/*_merged.fastq.gz')
                    .map { merged -> 
                        def sample = merged.simpleName.replaceAll('_merged', '')
                        def profile = file("${projectDir}/testdata/precomputed_metaphlan/metaphlan/\${sample}_profile.tsv")
                        tuple(sample, merged, profile)
                    }
                    .first()
                input[1] = file('${projectDir}/testdata/humann_demo_db')
                input[2] = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
                input[3] = '${projectDir}/testdata/metaphlan_db'
                """
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("HUMANN_RUN - generate functional profiles (Docker)") {
        config "conf/test.docker.config"
        
        when {
            params {
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/*_merged.fastq.gz')
                    .map { merged -> 
                        def sample = merged.simpleName.replaceAll('_merged', '')
                        def profile = file("${projectDir}/testdata/precomputed_metaphlan/metaphlan/\${sample}_profile.tsv")
                        tuple(sample, merged, profile)
                    }
                    .first()
                input[1] = file('${projectDir}/testdata/humann_demo_db')
                input[2] = 'mpa_vJan25_CHOCOPhlAnSGB_202503'
                input[3] = '${projectDir}/testdata/metaphlan_db'
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test HUMANN_INIT Process"
    script "tests/modules/humann_test_workflow.nf"
    workflow "test_humann_init"
    
    test("HUMANN_INIT - symlink existing database") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                // No input needed
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test HUMANN_MERGE Process"
    script "tests/modules/humann_test_workflow.nf"
    workflow "test_humann_merge"
    
    test("HUMANN_MERGE - merge outputs from multiple samples") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                input[0] = Channel.fromPath('${projectDir}/testdata/humann_outputs/*_genefamilies.tsv')
                    .collect()
                input[1] = Channel.fromPath('${projectDir}/testdata/humann_outputs/*_pathabundance.tsv')
                    .collect()
                input[2] = Channel.fromPath('${projectDir}/testdata/humann_outputs/*_pathcoverage.tsv')
                    .collect()
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test HUMANN_RENORM Process"
    script "tests/modules/humann_test_workflow.nf"
    workflow "test_humann_renorm"
    
    test("HUMANN_RENORM - renormalize to CPM") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                input[0] = file('${projectDir}/testdata/humann_outputs/humann_genefamilies_rpk.tsv')
                input[1] = file('${projectDir}/testdata/humann_outputs/humann_pathabundance_rpk.tsv')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}

nextflow_workflow {

    name "Test HUMANN_REGROUP Process"
    script "tests/modules/humann_test_workflow.nf"
    workflow "test_humann_regroup"
    
    test("HUMANN_REGROUP - regroup gene families") {
        config "conf/test.apptainer.config"
        
        when {
            params {
                humann_image = 'docker://biobakery/humann:latest'
            }
            workflow {
                """
                input[0] = file('${projectDir}/testdata/humann_outputs/humann_genefamilies_cpm.tsv')
                input[1] = file('${projectDir}/testdata/humann_demo_db')
                """
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

echo ""
echo "✓ Module test infrastructure updated!"
echo ""
echo "Test the changes with:"
echo "  java -jar ~/.nf-test/nf-test.jar test tests/modules/mergereads.nf.test"
echo "  java -jar ~/.nf-test/nf-test.jar test tests/modules/kneaddata.nf.test"
echo "  java -jar ~/.nf-test/nf-test.jar test tests/modules/metaphlan.nf.test"
echo "  java -jar ~/.nf-test/nf-test.jar test tests/modules/humann.nf.test"
echo ""
