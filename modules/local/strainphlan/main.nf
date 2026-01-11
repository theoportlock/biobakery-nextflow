process STRAINPHLAN_EXTRACT_MARKERS {
    label "process_medium"

    tag "extract markers from $sample"
    publishDir "$params.outdir/strainphlan/markers", pattern: "*.pkl", mode: "copy"
    container "$params.metaphlan_image"

    input:
    tuple val(sample), path(samout)

    output:
    tuple val(sample), path("${sample}.pkl"), emit: markers

    script:
    """
    sample2markers.py \
        --input ${samout} \
        --output_dir . \
        --nprocs ${task.cpus}
    """
}

process STRAINPHLAN_EXTRACT_REF_MARKERS {
    label "process_medium"

    tag "extract markers from reference genomes"
    publishDir "$params.outdir/strainphlan/reference_markers", pattern: "*.fna", mode: "copy"
    container "$params.metaphlan_image"

    input:
    path reference_genomes

    output:
    path "*.fna", emit: ref_markers

    script:
    """
    for genome in ${reference_genomes}/*.fna; do
        extract_markers.py \
            --input \$genome \
            --output_dir . \
            --nprocs ${task.cpus}
    done
    """
}

process STRAINPHLAN_BUILD_TREES {
    label "process_high"
    label "process_long"

    tag "build tree for $sgb"
    publishDir "$params.outdir/strainphlan/${sgb}", pattern: "*", mode: "copy"
    container "$params.metaphlan_image"

    input:
    val sgb
    path markers
    path ref_markers, stageAs: 'ref_markers/*'

    output:
    tuple val(sgb), path("${sgb}/RAxML_bestTree.${sgb}.StrainPhlAn4.tre"), emit: tree
    tuple val(sgb), path("${sgb}/*.aln"), emit: alignments
    tuple val(sgb), path("${sgb}/*.info"), emit: info
    path "${sgb}/*"

    script:
    def ref_opt = ref_markers.name != 'NO_REF' ? "--references ref_markers/*.fna" : ""
    """
    mkdir -p strainphlan_output
    
    strainphlan \
        --samples ${markers}/*.pkl \
        --output_dir strainphlan_output \
        --clades ${sgb} \
        --nprocs ${task.cpus} \
        --marker_in_n_samples ${params.min_samples_per_marker} \
        --sample_with_n_markers ${params.min_markers_per_sample} \
        --phylophlan_mode accurate \
        --mutation_rates \
        ${ref_opt}
    
    # Move output to expected location
    mv strainphlan_output/${sgb} .
    """
}

process STRAINPHLAN_TRANSMISSION {
    label "process_single"

    tag "transmission analysis for $sgb"
    publishDir "$params.outdir/strainphlan/${sgb}", pattern: "*.info", mode: "copy"
    publishDir "$params.outdir/strainphlan/${sgb}", pattern: "*.dist", mode: "copy"
    container "$params.metaphlan_image"

    input:
    tuple val(sgb), path(tree)
    path metadata

    output:
    tuple val(sgb), path("transmission_events.info"), emit: transmission
    path "${tree.baseName}.dist", optional: true, emit: distances

    script:
    def save_dist_opt = params.save_strain_distances ? "--save_dist" : ""
    """
    strain_transmission.py \
        --tree ${tree} \
        --metadata ${metadata} \
        --output_dir . \
        --threshold ${params.strain_threshold} \
        ${save_dist_opt}
    """
}

process STRAINPHLAN_AGGREGATE_RESULTS {
    label "process_single"

    tag "aggregate transmission events"
    publishDir "$params.outdir/strainphlan", mode: "copy"

    input:
    path transmission_files

    output:
    path "strainphlan_transmission_summary.csv"

    script:
    """
    #!/usr/bin/env python3
    import os
    import csv
    
    # Parse all transmission_events.info files
    results = []
    
    for info_file in "${transmission_files}".split():
        if not info_file.endswith('transmission_events.info'):
            continue
            
        # Extract SGB from parent directory name
        sgb = os.path.basename(os.path.dirname(info_file))
        
        with open(info_file, 'r') as f:
            lines = f.readlines()
            
            # Parse threshold
            threshold = None
            if len(lines) > 0 and 'threshold:' in lines[0]:
                threshold = lines[0].split(':')[1].strip()
            
            # Parse number of events
            n_events = None
            if len(lines) > 1 and 'Number of transmission events:' in lines[1]:
                n_events = lines[1].split(':')[1].strip()
            
            # Parse transmission pairs (skip header lines)
            pairs = []
            for line in lines[3:]:
                line = line.strip()
                if '<->' in line:
                    pair = line.split(' <-> ')
                    if len(pair) == 2:
                        pairs.append((pair[0], pair[1]))
            
            # Add to results
            for pair in pairs:
                results.append({
                    'SGB': sgb,
                    'Sample1': pair[0],
                    'Sample2': pair[1],
                    'Threshold': threshold,
                    'Status': 'transmission_event'
                })
            
            # If no pairs, still record the SGB
            if len(pairs) == 0:
                results.append({
                    'SGB': sgb,
                    'Sample1': 'NA',
                    'Sample2': 'NA',
                    'Threshold': threshold,
                    'Status': 'no_transmission'
                })
    
    # Write summary CSV
    with open('strainphlan_transmission_summary.csv', 'w', newline='') as csvfile:
        if len(results) > 0:
            writer = csv.DictWriter(csvfile, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        else:
            # Empty file with header only
            writer = csv.DictWriter(csvfile, fieldnames=['SGB', 'Sample1', 'Sample2', 'Threshold', 'Status'])
            writer.writeheader()
    """
}
