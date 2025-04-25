// General parameters
params.barcodes = "barcodes.txt"
params.fastq_dir = ""
params.regions_bed = "${workflow.projectDir}/ReferenceSequences/SIVregions_wBarcode.bed"
params.reference = "${workflow.projectDir}/ReferenceSequences/SIVMac239FullGenome_wBarcode.fas"
params.outdir = "${workflow.projectDir}/nf-results"
params.virus = "SIV"
params.min_read_length = -1
params.max_read_length = -1
params.min_depth = 1000
params.count_barcodes = false
params.gpu = false

process concatenateFastq {
    
    executor 'local'
    // Conda not needed, but initialized now so that the next process doesn't time out while creating environment
    conda "${workflow.projectDir}/medaka.yaml"

    input:
    path fastq_dir
    tuple val(barcode), val(sample_id), val(fragments)
    val outdir
    errorStrategy 'retry'
    maxRetries 2

    output:
    path("${sample_id}_concatenated.fastq.gz")

    shell:
    """
    echo "Sample ${sample_id}:"
    # Create any missing directories
    if [ ! -d ${outdir} ]; then
        mkdir ${outdir};
    fi
    if [ ! -d ${outdir}/${sample_id} ]; then
        mkdir ${outdir}/${sample_id};
    elif [ -d ${outdir}/${sample_id} ]; then
        rm -rf ${outdir}/${sample_id};
        mkdir ${outdir}/${sample_id};
    fi

    # Either concatenate all fastqs in a folder or just take the one
    if [ -d ${fastq_dir}/${barcode} ]; then
        echo "Concatenating fastq files..."
        cat ${fastq_dir}/${barcode}/*.fastq.gz > ${sample_id}_concatenated.fastq.gz;
    elif [ -f ${fastq_dir}/${sample_id}.fastq.gz ]; then
        cat ${fastq_dir}/${sample_id}.fastq.gz > ${sample_id}_concatenated.fastq.gz;
    fi
    cp ${sample_id}_concatenated.fastq.gz ${outdir}/${sample_id}/
    if [ ! -f "${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz" ]; then
        exit 1
    fi

    """

}

process qualityFilter {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 1
    time = { 120.minute * task.attempt}
    memory = { 4.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("${sample_id}_concatenated.fastq.gz")
    val outdir
    val min_read_length
    val max_read_length

    output:
    path("${sample_id}_filtered.fastq.gz")

    shell:
    """
    # Split processing by fragment
    echo "Sample ${sample_id}:"
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        # Create a file for each fragment
        mkdir -p "${outdir}/${sample_id}/\$i"
        length=\$(awk -v pat=\$i '\$1 == pat {print \$2}' "${workflow.projectDir}/ReferenceSequences/SIV_frag_sizes.txt")

        if [ ${min_read_length} -gt -1 ]; then
            min_read=${min_read_length}
        else
            min_read=\$((\${length} - 100))
        fi
        if [ ${max_read_length} -gt -1 ]; then
            max_read=${max_read_length}
        else
            max_read=\$((\${length} + 100))
        fi

        fastplong \
            -i "${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz" \
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
        mkdir -p "${outdir}/${sample_id}/\$i/fastplong"
        cp "${sample_id}_filtered.fastq.gz" "${outdir}/${sample_id}/\$i/fastplong/"
        cp "${sample_id}_fastp_report.html" "${outdir}/${sample_id}/\$i/fastplong/"
        cp "${sample_id}_fastp_report.json" "${outdir}/${sample_id}/\$i/fastplong/"
    done

    """
}

process flyeAssembly {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 2
    time = { 20.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("${sample_id}_filtered.fastq.gz")
    path reference
    val outdir

    output:
    path("flye_out/assembly.fasta")
    
    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        length=\$(awk -v pat=\$i '\$1 == pat {print \$2}' "${workflow.projectDir}/ReferenceSequences/SIV_frag_sizes.txt")

        { # Try to use assembly
            flye \
                --nano-raw "${outdir}/${sample_id}/\$i/fastplong/${sample_id}_filtered.fastq.gz" \
                --out-dir flye_out \
                --genome-size \${length} \
                --min-overlap 1000 \
                --asm-coverage 50 \
                --iterations 5 \
                --no-alt-contigs \
                --threads 8
        } || { # Use reference if assebmly can't be built
                mkdir -p flye_out
                cat ${reference} > flye_out/assembly.fasta
        }
        cp -r flye_out/ ${outdir}/${sample_id}/\$i
    done
    """

}

process cleanContigs {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 1
    time = { 5.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("flye_out/assembly.fasta")
    path reference
    path regions_bed
    val outdir

    output:
    path("${sample_id}_assembly.fasta")
    
    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})

    # Split reference into separate contigs
    seqtk subseq ${reference} ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/scripts/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}
    

    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        echo "\$i" > frags.txt
        seqtk subseq renamed_contigs.fasta frags.txt > ref.fasta
        python "${workflow.projectDir}/scripts/reorder_contigs.py" \
            -r ref.fasta \
            -a "${outdir}/${sample_id}/\$i/flye_out/assembly.fasta" \
            -o "${sample_id}_assembly.fasta" \
            --concat
        cp "${sample_id}_assembly.fasta" "${outdir}/${sample_id}/\$i/"
    done
    """
}

process polishAssembly {

    executor 'slurm'
    if(params.gpu) {
        conda "${workflow.projectDir}/medaka.yaml"
        clusterOptions = '-A b1042 --gres=gpu:a100:1'
        queue = 'genomics-gpu'
    }
    else {
        conda "${workflow.projectDir}/medaka.yaml"
        clusterOptions = '-A b1042'
        queue = 'genomics'
    }
    cpus = 1
    time = { 30.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("${sample_id}_assembly.fasta")
    path("${sample_id}_filtered.fastq.gz")
    val outdir
    val gpu

    output:
    path("medaka_out/consensus.fasta")

    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        medaka_consensus \
            -i "${outdir}/${sample_id}/\$i/fastplong/${sample_id}_filtered.fastq.gz" \
            -d "${outdir}/${sample_id}/\$i/${sample_id}_assembly.fasta" \
            -m r1041_e82_400bps_sup_v5.0.0 \
            -o medaka_out \
            -t 8
        cp medaka_out/consensus.fasta "${outdir}/${sample_id}/\$i/${sample_id}_consensus.fasta"
    done
    """

}

process alignReads {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 2
    time = { 10.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("${sample_id}_filtered.fastq.gz")
    path("medaka_out/consensus.fasta")
    val outdir

    output:
    tuple path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")

    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})
    
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        minimap2 -ax lr:hq "${outdir}/${sample_id}/\$i/${sample_id}_consensus.fasta" "${outdir}/${sample_id}/\$i/fastplong/${sample_id}_filtered.fastq.gz" | \
            samtools sort -o "${sample_id}_aligned.bam"
        samtools index "${sample_id}_aligned.bam"
        cp "${sample_id}_aligned.bam" "${outdir}/${sample_id}/\$i/"
        cp "${sample_id}_aligned.bam.bai" "${outdir}/${sample_id}/\$i/"
    done
    """

}

process trimConsensus {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 1
    time = { 5.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("medaka_out/consensus.fasta")
    tuple path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")
    val virus
    val outdir
    val min_depth

    output:
    path("${sample_id}_consensus_trimmed.fasta")
    
    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})

    for i in \${frags[@]}; do

        echo "Trimming \${i}..."
        python "${workflow.projectDir}/scripts/trim_contigs_by_depth.py" \
            -f "${outdir}/${sample_id}/\$i/${sample_id}_consensus.fasta" \
            -b "${outdir}/${sample_id}/\$i/${sample_id}_aligned.bam" \
            -o "${sample_id}_consensus_trimmed.fasta" \
            -d ${min_depth}

        if [ ${virus} == "SIV" ] && [ \$i == "Fragment_3" ]; then
            python "${workflow.projectDir}/scripts/split_fasta_by_delimiter.py" \
                -i "${sample_id}_consensus_trimmed.fasta" \
                -d "acgcgagccatgnnnnnnnnnncgatgcgcgcgt" \
                -o "${sample_id}_consensus_split.fasta"
            rm "${sample_id}_consensus_trimmed.fasta"
            mv "${sample_id}_consensus_split.fasta" "${sample_id}_consensus_trimmed.fasta"
        fi

        cp "${sample_id}_consensus_trimmed.fasta" "${outdir}/${sample_id}/\$i/"

        echo "Realigning \${i}..."
        minimap2 -ax lr:hq "${sample_id}_consensus_trimmed.fasta" "${outdir}/${sample_id}/\$i/fastplong/${sample_id}_filtered.fastq.gz" | \
            samtools sort -o "${sample_id}_aligned.bam"
        samtools index "${sample_id}_aligned.bam"
        cp "${sample_id}_aligned.bam" "${outdir}/${sample_id}/\$i/"
        cp "${sample_id}_aligned.bam.bai" "${outdir}/${sample_id}/\$i/"

    done
    """

}

process haplotypes {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 2
    time = { 10.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    tuple path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai")
    path reference
    path regions_bed
    val outdir
    val virus

    output:
    path("read_counts.txt")

    shell:
    """
    echo "Sample ${sample_id}:"
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        
        # Count each 'haplotype'
        python "${workflow.projectDir}/scripts/contig_read_counter.py" "${outdir}/${sample_id}/\$i/${sample_id}_aligned.bam" --output-bam "${sample_id}_filtered.bam"

        cp "read_counts.txt" "${outdir}/${sample_id}/\$i/"
    done
    """

}

process countBarcodes {

    executor 'slurm'
    conda "${workflow.projectDir}/medaka.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 2
    time = { 8.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path("${sample_id}_concatenated.fastq.gz")
    path outdir
    val gpu

    shell:
    """
    if [ ! -d ${outdir}/${sample_id} ]; then
        mkdir ${outdir}/${sample_id};
    fi
    if [ ! -d ${outdir}/${sample_id}/barcode/ ]; then
        mkdir ${outdir}/${sample_id}/barcode/
    fi

    minimap2 -ax lr:hq ${workflow.projectDir}/ReferenceSequences/barcode_reference.fas ${sample_id}_concatenated.fastq.gz > ${sample_id}_barcode.sam
    samtools view -F 4 -bS -h -O BAM -e 'rlen>233' ${sample_id}_barcode.sam > ${sample_id}_barcode.bam
    samtools sort ${sample_id}_barcode.bam > ${sample_id}_barcode_sorted.bam; samtools index ${sample_id}_barcode_sorted.bam
    samtools depth -a -r "barcode" -@ 4 ${sample_id}_barcode_sorted.bam > ${outdir}/${sample_id}/barcode/barcode_depth.txt
    perl ${workflow.projectDir}/scripts/bam_barcode_count.pl -b ${sample_id}_barcode_sorted.bam -s 113 -e 122 > ${outdir}/${sample_id}/barcode/barcodes.txt
    cp ${sample_id}_barcode_sorted.bam ${outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam
    cp ${sample_id}_barcode_sorted.bam.bai ${outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam.bai

    """
}

workflow {
    
    // Read in csv
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv(header:false)
        .map { row -> tuple(row[0], row[1], row[2]) }
        .set { samples_ch }
    
    // What fragments contain a barcode
    def barcoded_regions = ["Fragment_3", "FL"]

    // Filter for only samples with a barcode - will be used in countBarcodes process
    samples_ch
        .filter { barcode, sample_id, fragments -> fragments != null && barcoded_regions.any { fragments.contains(it) }}
        .set { barcoded_ch }
    
    // Run workflow
    concatenated = concatenateFastq(params.fastq_dir, samples_ch, params.outdir)
    quality = qualityFilter(samples_ch, concatenated, params.outdir, params.min_read_length, params.max_read_length)
    flye = flyeAssembly(samples_ch, quality, params.reference, params.outdir)
    cleaned = cleanContigs(samples_ch, flye, params.reference, params.regions_bed, params.outdir)
    polished = polishAssembly(samples_ch, cleaned, quality, params.outdir, params.gpu)
    aligned = alignReads(samples_ch, quality, polished, params.outdir)
    trimmed = trimConsensus(samples_ch, polished, aligned, params.virus, params.outdir, params.min_depth)
    haplotypes(samples_ch, aligned, params.reference, params.regions_bed, params.outdir, params.virus)
    // Optional process to split/count barcoded fragments
    if (params.count_barcodes)
        countBarcodes(barcoded_ch, concatenated, params.outdir, params.gpu)
}
