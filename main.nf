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
        echo "Concatenating file failed"
        exit 1
    fi

    """

}

// 2. Quality filter
process qualityFilter {

    input:
    tuple val(barcode), val(sample_id), path(concat_fastq)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_filtered.fastq.gz")

    script:
    """

    echo "Sample ${sample_id}"

    # Create a file sample
    mkdir -p "${params.outdir}/${sample_id}"

    if [ ${params.min_read_length} -gt -1 ]; then
        min_read=${params.min_read_length}
    else
        min_read=1200
    fi
    if [ ${params.max_read_length} -gt -1 ]; then
        max_read=${params.max_read_length}
    else
        max_read=20000
    fi

    fastplong \
        -i "${sample_id}_concatenated.fastq.gz" \
        -o "${sample_id}_filtered.fastq.gz" \
        --qualified_quality_phred 15 \
        --length_required \$min_read \
        --length_limit \$max_read \
        --thread 4 \
        --html "${sample_id}_fastp_report.html" \
        --json "${sample_id}_fastp_report.json"
        
    echo "Finished filtering"
    echo "Moving files..."
    mkdir -p "${params.outdir}/${sample_id}/fastplong"
    cp -f "${sample_id}_filtered.fastq.gz" "${params.outdir}/${sample_id}/fastplong/"
    cp -f "${sample_id}_fastp_report.html" "${params.outdir}/${sample_id}/fastplong/"
    cp -f "${sample_id}_fastp_report.json" "${params.outdir}/${sample_id}/fastplong/"

    """
}

// ----------------------------------------------
// -------------- De Novo -----------------------
// ----------------------------------------------
// 3. Flye assembly
process flyeAssembly {

    input:
    tuple val(barcode), val(sample_id), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), path("flye_out/assembly.fasta")
    
    script:
    """
    echo "Sample ${sample_id}:"

    # Try to use assembly
    run_flye () {
        flye \
            --nano-raw "${sample_id}_filtered.fastq.gz" \
            --out-dir flye_out \
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

    cp -rf flye_out/ "${params.outdir}/${sample_id}/"

    """

}

//4. Polish Assembly
process polishAssembly {

    input:
    tuple val(barcode), val(sample_id), path(assembly_fasta), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), path("medaka_out/consensus.fasta")

    script:
    """
    echo "Sample ${sample_id}:"

    medaka_consensus \
        -i "${filtered_fastq}" \
        -d "${assembly_fasta}" \
        -m r1041_e82_400bps_sup_v5.0.0 \
        -o medaka_out \
        -t 8 \
        -f -x
    cp -f medaka_out/consensus.fasta "${params.outdir}/${sample_id}/${sample_id}_consensus.fasta"

    """

}

// 5. Align reads back
process alignReads {

    input:
    tuple val(barcode), val(sample_id), path(consensus_fasta), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")

    script:
    """
    echo "Sample ${sample_id}:"
    
    minimap2 -ax lr:hq "consensus.fasta" "${sample_id}_filtered.fastq.gz" | \
        samtools sort -o "${sample_id}_aligned.bam"
    samtools index "${sample_id}_aligned.bam"
    cp -f "${sample_id}_aligned.bam" "${params.outdir}/${sample_id}/"
    cp -f "${sample_id}_aligned.bam.bai" "${params.outdir}/${sample_id}/"

    """

}

// ---------------------------------------------------
// ---------------- Reference-based ------------------
// ---------------------------------------------------

// 3. First Alignment
process firstAlignment {
    
    input:
    tuple val(barcode), val(sample_id), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_aligned.bam")

    script:
    """
    echo "Sample ${sample_id}:"

    minimap2 -ax lr:hq "${params.reference}" "${sample_id}_filtered.fastq.gz" > aligned.sam
    samtools view -bS -F 0x904 -h aligned.sam > aligned.bam
    samtools sort aligned.bam > "${sample_id}_aligned.bam"
    """
}

// 4. First Consensus
process firstConsensus {

    input:
    tuple val(barcode), val(sample_id), path(aligned_bam)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_con1.fa")

    script:
    """
    samtools mpileup -aa -A -Q 0 ${aligned_bam} | ivar consensus -p "${sample_id}_con1" -q 10 -m 100 -k
    """
}

// 5. Realign
process realign {

    input:
    tuple val(barcode), val(sample_id), path(firstCon), path(filtered_fastq)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}.bam")

    script:
    """
    minimap2 -ax lr:hq "${firstCon}" "${filtered_fastq}" > aligned.sam
    samtools view -bS -F 4 -h aligned.sam > aligned.bam
    samtools sort aligned.bam > "${sample_id}.bam"
    samtools index "${sample_id}.bam"

    cp -f "${sample_id}.bam" "${params.outdir}/${sample_id}/"
    cp -f "${sample_id}.bam.bai" "${params.outdir}/${sample_id}/"
    """
}

// 6. Final consensus
process finalConsensus {

    input:
    tuple val(barcode), val(sample_id), path(aligned_bam)

    output:
    tuple val(barcode), val(sample_id), path("${sample_id}_consensus.fa")

    script:
    """
    samtools mpileup -aa -A -Q 0 "${aligned_bam}" | ivar consensus -p "${sample_id}_consensus" -q 10 -t 0.5 -n N

    cp -f "${sample_id}"_consensus.* "${params.outdir}/${sample_id}/"
    """
}

//-------------------------------------------------
// ------------------- Workflow -------------------
// ------------------------------------------------

workflow {
    
    // Read in csv
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv(header:false)
        .map { row -> tuple(row[0], row[1]) }
        .set { samples_ch }
    
    // Concatenate fastq once per sample
    concatenated = concatenateFastq(samples_ch)
    // Quality filter
    quality = qualityFilter(concatenated)

    // Force denovo if no reference given
    if ( params.reference == null || params.reference == "" ) {
        params.denovo = true
    }

    // Denovo Assembly
    if ( params.denovo ) {
        assembled = flyeAssembly(quality)

        // Join cleaned and quality for polishing
        polish_input = assembled
            .join(quality, by: [1, 1])
            .map { sample, barcode, assembly, _barcode1, filtered -> 
                tuple(barcode, sample, assembly, filtered)}
        polished = polishAssembly(polish_input)

        // Join polished and quality for alignment
        align_input = polished
            .join(quality, by: [1, 1])
            .map { sample, barcode, consensus, _barcode1, filtered -> 
                tuple(barcode, sample, consensus, filtered)}
        aligned = alignReads(align_input)
    }
    // Reference-guided Assembly
    else {
        align1 = firstAlignment(quality)
        con1 = firstConsensus(align1)

        realign_input = con1
            .join(quality, by: [1, 1])
            .map { sample, barcode, consensus, _barcode1, filtered ->
                tuple(barcode, sample, consensus, filtered)}
        realigned = realign(realign_input)
        conf = finalConsensus(realigned)
    }

}
