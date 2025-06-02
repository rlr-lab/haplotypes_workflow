#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// 1. Concatenate FASTQ
process concatenateFastq {
    
    input:
    tuple val(barcode), val(sample_id)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_concatenated.fastq.gz")

    script:
    """
    echo "Sample ${sample_id}:"
    # Create any missing directories
    if [ ! -d "${params.outdir}" ]; then
        mkdir "${params.outdir}"
    fi
    if [ ! -d "${params.outdir}/${sample_id}" ]; then
        mkdir "${params.outdir}/${sample_id}"
    elif [ -d "${params.outdir}/${sample_id}" ]; then
        rm -rf "${params.outdir}/${sample_id}"
        mkdir "${params.outdir}/${sample_id}"
    fi

    # Either concatenate all fastqs in a folder or just take the one
    if [ -d "${params.fastq_dir}/${barcode}" ]; then
        echo "Concatenating fastq files..."
        cat ${params.fastq_dir}/${barcode}/*.fastq.gz > "${sample_id}_concatenated.fastq.gz"
    elif [ -f "${params.fastq_dir}/${sample_id}.fastq.gz" ]; then
        cat "${params.fastq_dir}/${sample_id}.fastq.gz" > "${sample_id}_concatenated.fastq.gz"
    fi

    cp -f "${sample_id}_concatenated.fastq.gz" "${params.outdir}/${sample_id}/"

    if [ ! -f "${params.outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz" ]; then
        exit 1
    fi

    """

}

// 2. Quality filter
process qualityFilter {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(concat_fastq)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("${sample_id}_filtered.fastq.gz")

    script:
    """
    # Split processing by fragment
    echo "Sample ${sample_id}, ${fragments}:"

    # Create a file for each fragment
    mkdir -p "${params.outdir}/${sample_id}/${fragments}"
    length=\$(awk -v pat=${fragments} '\$4 == pat {print \$3 - \$2}' "${params.regions_bed}")

    if [ ${params.min_read_length} -gt -1 ]; then
        min_read=${params.min_read_length}
    else
        min_read=\$((\${length} - 100))
    fi
    if [ ${params.max_read_length} -gt -1 ]; then
        max_read=${params.max_read_length}
    else
        max_read=\$((\${length} + 100))
    fi

    fastplong \
        -i "${sample_id}_concatenated.fastq.gz" \
        -o "${sample_id}_filtered.fastq.gz" \
        --qualified_quality_phred 15 \
        --length_required \$min_read \
        --length_limit \$max_read \
        --adapter_fasta "${workflow.projectDir}/ReferenceSequences/SIVprimers.fasta" \
        --thread 4 \
        --html "${sample_id}_fastp_report.html" \
        --json "${sample_id}_fastp_report.json"
        
    echo "Finished filtering"
    echo "Moving files..."
    mkdir -p "${params.outdir}/${sample_id}/${fragments}/fastplong"
    cp -f "${sample_id}_filtered.fastq.gz" "${params.outdir}/${sample_id}/${fragments}/fastplong/"
    cp -f "${sample_id}_fastp_report.html" "${params.outdir}/${sample_id}/${fragments}/fastplong/"
    cp -f "${sample_id}_fastp_report.json" "${params.outdir}/${sample_id}/${fragments}/fastplong/"

    """
}

// 3. Flye assembly
process flyeAssembly {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("flye_out/assembly.fasta")
    
    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"

    length=\$(awk -v pat=${fragments} '\$4 == pat {print \$3 - \$2}' "${params.regions_bed}")

    # Try to use assembly
    run_flye () {
        flye \
            --nano-raw "${sample_id}_filtered.fastq.gz" \
            --out-dir flye_out \
            --genome-size \${length} \
            --min-overlap 1000 \
            --asm-coverage 50 \
            --iterations 5 \
            --no-alt-contigs \
            --threads 8
    }

    if run_flye; then
        echo "Flye consensus generated"
    else
        echo "Flye failed - copying reference"
        mkdir -p flye_out
        cat ${params.reference} > flye_out/assembly.fasta
    fi

    cp -rf flye_out/ "${params.outdir}/${sample_id}/${fragments}"

    """

}

// 4. Clean and reorder contigs
process cleanContigs {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(assembly_fasta)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("${sample_id}_assembly.fasta")
    
    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"

    # Split reference into separate contigs
    seqtk subseq ${params.reference} ${params.regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/scripts/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${params.regions_bed}
    
    echo ${fragments} > frags.txt
    seqtk subseq renamed_contigs.fasta frags.txt > ref.fasta
    python "${workflow.projectDir}/scripts/reorder_contigs.py" \
        -r ref.fasta \
        -a "assembly.fasta" \
        -o "${sample_id}_assembly.fasta" \
        --concat
    cp -f "${sample_id}_assembly.fasta" "${params.outdir}/${sample_id}/${fragments}/"

    """
}

//5. Polish Assembly
process polishAssembly {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(assembly_fasta), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("medaka_out/consensus.fasta")

    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"

    medaka_consensus \
        -i "${sample_id}_filtered.fastq.gz" \
        -d "${sample_id}_assembly.fasta" \
        -m r1041_e82_400bps_sup_v5.0.0 \
        -o medaka_out \
        -t 8 \
        -f -x
    cp -f medaka_out/consensus.fasta "${params.outdir}/${sample_id}/${fragments}/${sample_id}_consensus.fasta"

    """

}

// 6. Align reads back
process alignReads {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(consensus_fasta), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")

    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"
    
    minimap2 -ax lr:hq "consensus.fasta" "${sample_id}_filtered.fastq.gz" | \
        samtools sort -o "${sample_id}_aligned.bam"
    samtools index "${sample_id}_aligned.bam"
    cp -f "${sample_id}_aligned.bam" "${params.outdir}/${sample_id}/${fragments}/"
    cp -f "${sample_id}_aligned.bam.bai" "${params.outdir}/${sample_id}/${fragments}/"

    """

}

// 7. Trim consensus by depth
process trimConsensus {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(consensus_fasta), path(aligned_bam), path(aligned_bai)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("${sample_id}_consensus_trimmed.fasta"), path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")
    
    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"

    # Drop shortest 10% of aligned reads
    python "${workflow.projectDir}/scripts/filter_by_aligned_length.py" \
        "${sample_id}_aligned.bam" \
        length_filtered.bam
    samtools sort length_filtered.bam > length_filtered_sorted.bam
    samtools index length_filtered_sorted.bam

    depth=\$(samtools depth length_filtered_sorted.bam | awk '{sum+=\$3} END {print (sum/NR)*0.9}')
    echo "Trimming to a minimum depth of \${depth}"

    # Trim contig to remove bases on ends with low coverage
    python "${workflow.projectDir}/scripts/trim_contigs_by_depth.py" \
        -f "consensus.fasta" \
        -b length_filtered_sorted.bam \
        -o "${sample_id}_consensus_trimmed.fasta" \
        -d \${depth}

    # Split SIV Fragment_3 into 2 sequences around barcode
    if [ ${params.virus} == "SIV" ] && [ ${fragments} == "Fragment_3" ]; then
        splitFasta() {
            python "${workflow.projectDir}/scripts/split_fasta_by_delimiter.py" \
                -i "${sample_id}_consensus_trimmed.fasta" \
                -d "acgcgagccatgnnnnnnnnnncgatgcgcgcgt" \
                -o "${sample_id}_consensus_split.fasta"
        }
        if splitFasta && [ -f "${sample_id}_consensus_split.fasta" ]; then
            echo "Contig split by delimiter"
            echo "Further processing done over two contigs"
            rm "${sample_id}_consensus_trimmed.fasta"
            mv "${sample_id}_consensus_split.fasta" "${sample_id}_consensus_trimmed.fasta"
        else
            echo "Delimiter could not be found"
            echo "Assembly kept as one contig"
        fi
    fi

    # Rename consensus sequence in fasta
    cp "${sample_id}_consensus_trimmed.fasta" temp_rename.fasta
    python "${workflow.projectDir}/scripts/fastaRegexRemove.py" \
        -i temp_rename.fasta \
        -r "assembled_sequence" \
        -n "${fragments}" \
        -o "${sample_id}_consensus_trimmed.fasta"

    cp -f "${sample_id}_consensus_trimmed.fasta" "${params.outdir}/${sample_id}/${fragments}/"

    echo "Realigning..."
    minimap2 -ax lr:hq "${sample_id}_consensus_trimmed.fasta" "${params.outdir}/${sample_id}/${fragments}/fastplong/${sample_id}_filtered.fastq.gz" > aligned.sam
    samtools view -bSh -F 4 aligned.sam > aligned.bam
    samtools sort -o "${sample_id}_aligned.bam" aligned.bam
    samtools index "${sample_id}_aligned.bam"
    cp -f "${sample_id}_aligned.bam" "${params.outdir}/${sample_id}/${fragments}/"
    cp -f "${sample_id}_aligned.bam.bai" "${params.outdir}/${sample_id}/${fragments}/"

    """

}

// 8. Haplotype counting
process haplotypes {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(trimmed_fasta), path(aligned_bam), path(aligned_bai)

    output:
    tuple val(barcode), val(sample_id), val(fragments), path("read_counts.txt")

    script:
    """
    echo "Sample ${sample_id}, ${fragments}:"
        
    # Count each 'haplotype'
    python "${workflow.projectDir}/scripts/contig_read_counter.py" "${sample_id}_aligned.bam" \
        --output-bam "${sample_id}_filtered.bam" \
        --min-mapping-quality 15

    cp -f "read_counts.txt" "${params.outdir}/${sample_id}/${fragments}/"

    """

}

process countBarcodes {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path("${sample_id}_filtered.fastq.gz")

    script:
    """

    if [ ! -d "${params.outdir}/${sample_id}/barcode/" ]; then
        mkdir "${params.outdir}/${sample_id}/barcode/"
    fi

    minimap2 -ax lr:hq "${workflow.projectDir}/ReferenceSequences/barcode_reference.fas" "${sample_id}_filtered.fastq.gz" > "${sample_id}_barcode.sam"
    samtools view -F 4 -bS -h -O BAM -e 'rlen>233' "${sample_id}_barcode.sam" > "${sample_id}_barcode.bam"
    samtools sort "${sample_id}_barcode.bam" > "${sample_id}_barcode_sorted.bam"; samtools index "${sample_id}_barcode_sorted.bam"
    samtools depth -a -r "barcode" -@ 4 "${sample_id}_barcode_sorted.bam" > "${params.outdir}/${sample_id}/barcode/barcode_depth.txt"
    perl "${workflow.projectDir}/scripts/bam_barcode_count.pl" -b "${sample_id}_barcode_sorted.bam" -s 113 -e 122 > "${params.outdir}/${sample_id}/barcode/barcodes.txt"
    cp -f "${sample_id}_barcode_sorted.bam" "${params.outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam"
    cp -f "${sample_id}_barcode_sorted.bam.bai" "${params.outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam.bai"

    """
}

process RVHaplo {

    input:
    tuple val(barcode), val(sample_id), val(fragments), path(trimmed_fasta), path(aligned_bam), path(aligned_bai)

    script:
    """
    # Run RVHaplo

    echo "Sample ${sample_id}, ${fragments}:"
 
    RVHaploPath=${workflow.projectDir}"/scripts/rvhaplo.sh"
 
    "\${RVHaploPath}" \
        -i "${sample_id}_aligned.bam" \
        -r "${sample_id}_consensus_trimmed.fasta" \
        -o ${params.outdir}/${sample_id}/${fragments} \
        -p ${sample_id} \
        -t 8 -e 0.1 -sg 10 -a 0.001 -ss 5
    
    echo "Done"
    """

}

workflow {
    
    // Read in csv
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv(header:false)
        .map { row -> tuple(row[0], row[1]) }
        .set { samples_ch }
    
    // Concatenate fastq once per sample
    concatenated = concatenateFastq(samples_ch)

    // Split to flatmap of fragments
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv(header:false)
        .flatMap { row ->
            def barcode = row[0]
            def sample_id = row[1]
            def fragments = row[2].tokenize(' ')
            fragments.collect { fragment -> tuple(barcode, sample_id, fragment) }
        }
        .set { fragmented_ch }
    
    // Combine with concatenated file
    quality_input = fragmented_ch
        .combine(concatenated, by: 1)
        .map { sample, barcode, fragment, _barcode1, concat ->
            tuple(barcode, sample, fragment, concat)
        }

    // Continue with workflow
    quality = qualityFilter(quality_input)
    assembled = flyeAssembly(quality)
    cleaned = cleanContigs(assembled)

    // Join cleaned and quality for polishing
    polish_input = cleaned
        .join(quality, by: [1, 1])
        .map { sample, barcode, fragments, assembly, _barcode1, _fragments1, filtered -> 
            tuple(barcode, sample, fragments, assembly, filtered)}
    polished = polishAssembly(polish_input)

    // Join polished and quality for alignment
    align_input = polished
        .join(quality, by: [1, 1])
        .map { sample, barcode, fragments, consensus, _barcode1, _fragments1, filtered -> 
            tuple(barcode, sample, fragments, consensus, filtered)}
    aligned = alignReads(align_input)

    // Join polished and aligned for trimming
    trim_input = polished
        .join(aligned, by: [1, 1])
        .map { sample, barcode, fragments, consensus, _barcode1, _fragments1, bam, bai -> 
            tuple(barcode, sample, fragments, consensus, bam, bai)}
    trimmed = trimConsensus(trim_input)

    // Haplotyping
    if ( params.rvhaplo ) {
        RVHaplo(trimmed)
    }
    else {
        haplotypes(trimmed)
    }

    // Optional process to split/count barcoded fragments
    if (params.count_barcodes) {
        // What fragments contain a barcode
        def barcoded_regions = ["Fragment_3", "FL"]

        // Filter input for barcoded samples, join with quality output
        barcoded_ch = quality
            .filter { _barcode, _sample_id, fragments, _path -> fragments != null && barcoded_regions.any { fragments.contains(it) }}
        
        countBarcodes(barcoded_ch)
    }
}
