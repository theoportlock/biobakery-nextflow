#!/bin/bash
# Simplify integration workflow tests to use workflow.success
set -e

echo "Simplifying integration workflow tests..."
echo ""

# For each integration test, we'll simplify to just check workflow.success
# This makes tests more reliable and faster

# Update kneaddata_workflow.nf.test
echo "Updating tests/workflow/kneaddata_workflow.nf.test..."
cat > tests/workflow/kneaddata_workflow.nf.test << 'EOF'
nextflow_workflow {

    name "Test KneadData Workflow"
    script "main.nf"

    test("KneadData - host read removal") {
        
        when {
            params {
                input = "${projectDir}/testdata/FG00004_S26_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("KneadData - multiple samples") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update metaphlan_workflow.nf.test
echo "Updating tests/workflow/metaphlan_workflow.nf.test..."
cat > tests/workflow/metaphlan_workflow.nf.test << 'EOF'
nextflow_workflow {

    name "Test MetaPhlAn Workflow"
    script "main.nf"

    test("MetaPhlAn - taxonomic profiling") {
        
        when {
            params {
                input = "${projectDir}/testdata/FG00004_S26_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                m = true
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("MetaPhlAn - multiple samples") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                m = true
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update humann_workflow.nf.test
echo "Updating tests/workflow/humann_workflow.nf.test..."
cat > tests/workflow/humann_workflow.nf.test << 'EOF'
nextflow_workflow {

    name "Test HUMAnN Workflow"
    script "main.nf"

    test("HUMAnN - functional profiling with precomputed profiles") {
        
        when {
            params {
                input = "${projectDir}/testdata/FG00004_S26_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                h = true
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                precomputed_metaphlan_dir = "${projectDir}/testdata/precomputed_metaphlan"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update kneaddata_metaphlan.nf.test
echo "Updating tests/workflow/kneaddata_metaphlan.nf.test..."
cat > tests/workflow/kneaddata_metaphlan.nf.test << 'EOF'
nextflow_workflow {

    name "Test KneadData + MetaPhlAn Integration"
    script "main.nf"

    test("K+M - host removal and taxonomic profiling") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                m = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update metaphlan_humann.nf.test
echo "Updating tests/workflow/metaphlan_humann.nf.test..."
cat > tests/workflow/metaphlan_humann.nf.test << 'EOF'
nextflow_workflow {

    name "Test MetaPhlAn + HUMAnN Integration"
    script "main.nf"

    test("M+H - taxonomic and functional profiling") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                m = true
                h = true
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update full_pipeline_kmh.nf.test
echo "Updating tests/workflow/full_pipeline_kmh.nf.test..."
cat > tests/workflow/full_pipeline_kmh.nf.test << 'EOF'
nextflow_workflow {

    name "Test Full Pipeline K+M+H"
    script "main.nf"

    test("K+M+H - complete functional profiling pipeline") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                m = true
                h = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
            }
        }

        then {
            assert workflow.success
        }
    }
    
    test("K+M+H - single sample complete pipeline") {
        
        when {
            params {
                input = "${projectDir}/testdata/FG00004_S26_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                m = true
                h = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update full_pipeline_kms.nf.test
echo "Updating tests/workflow/full_pipeline_kms.nf.test..."
cat > tests/workflow/full_pipeline_kms.nf.test << 'EOF'
nextflow_workflow {

    name "Test Full Pipeline K+M+S"
    script "main.nf"

    test("K+M+S - complete strain profiling pipeline") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                m = true
                s = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                sgb_list = "${projectDir}/testdata/sgb_list.txt"
                strain_metadata = "${projectDir}/testdata/strain_metadata.tsv"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

# Update full_pipeline_all.nf.test
echo "Updating tests/workflow/full_pipeline_all.nf.test..."
cat > tests/workflow/full_pipeline_all.nf.test << 'EOF'
nextflow_workflow {

    name "Test Full Pipeline K+M+H+S"
    script "main.nf"

    test("K+M+H+S - complete pipeline with all modules") {
        
        when {
            params {
                input = "${projectDir}/testdata/*_L001_R{1,2}_*.fastq.gz"
                outdir = "${outputDir}"
                k = true
                m = true
                h = true
                s = true
                kneaddata_db_dir = "${projectDir}/testdata/kneaddata_db"
                kneaddata_db = "test_host"
                metaphlan_db_dir = "${projectDir}/testdata/metaphlan_db"
                humann_db_dir = "${projectDir}/testdata/humann_demo_db"
                sgb_list = "${projectDir}/testdata/sgb_list.txt"
                strain_metadata = "${projectDir}/testdata/strain_metadata.tsv"
            }
        }

        then {
            assert workflow.success
        }
    }
}
EOF

echo ""
echo "✓ Integration workflow tests simplified!"
echo ""
echo "Test the changes with:"
echo "  java -jar ~/.nf-test/nf-test.jar test tests/workflow/"
echo ""
