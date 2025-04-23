// General parameters
params.barcodes = "barcodes.txt"
params.fastq_dir = ""
params.regions_bed = "${workflow.projectDir}/ReferenceSequences/SIVregions_wBarcode.bed"
params.reference = "${workflow.projectDir}/ReferenceSequences/SIVMac239FullGenome_wBarcode.fas"
params.outdir = "${workflow.projectDir}/nf-results"
params.virus = "SIV"
params.min_read_length = 1200
params.min_depth = 1000
params.split_barcode = false
params.gpu = false
// RVHaplo parameters
params.subgraphs = 1
params.abundance = 0.001
params.smallest_snv = 5

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
    path "${sample_id}_concatenated.fastq.gz"

    shell:
    """

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
    if [ ! -f ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz ]; then
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
    time = { 45.minute * task.attempt}
    memory = { 4.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path "${sample_id}_concatenated.fastq.gz"
    val outdir

    output:
    path "${sample_id}_filtered.fastq.gz"

    shell:
    """
    # Split processing by fragment
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        # Create a file for each fragment
        mkdir -p "${outdir}/${sample_id}/\$i"
        length=\$(awk -v pat=\$i '\$1 == pat {print \$2}' "${workflow.projectDir}/ReferenceSequences/SIV_frag_sizes.txt")
        fastplong \
            -i "${sample_id}_concatenated.fastq.gz" \
            -o "${sample_id}_filtered.fastq.gz" \
            --qualified_quality_phred 15 \
            --length_required \$((\${length} - 100)) \
            --length_limit \$((\${length} + 100)) \
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
    cpus = 4
    time = { 15.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path "${sample_id}_filtered.fastq.gz"
    val outdir

    output:
    path "flye_out/assembly.fasta"
    
    shell:
    """
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        length=\$(awk -v pat=\$i '\$1 == pat {print \$2}' "${workflow.projectDir}/ReferenceSequences/SIV_frag_sizes.txt")
        flye \
            --nano-raw "${outdir}/${sample_id}/\$i/fastplong/${sample_id}_filtered.fastq.gz" \
            --out-dir flye_out \
            --genome-size \${length} \
            --min-overlap 1000 \
            --asm-coverage 50 \
            --iterations 5 \
            --no-alt-contigs \
            --threads 8
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
    path "flye_out/assembly.fasta"
    path reference
    path regions_bed
    val outdir

    output:
    path "${sample_id}_assembly.fasta"
    
    shell:
    """
    frags=(${fragments})

    # Split reference into separate contigs
    seqtk subseq ${reference} ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}
    

    for i in \${frags[@]}; do
        echo "Processing \${i}..."
        echo "\$i" > frags.txt
        seqtk subseq renamed_contigs.fasta frags.txt > ref.fasta
        python "${workflow.projectDir}/reorder_contigs.py" \
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
    cpus = 4
    time = { 30.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path "${sample_id}_assembly.fasta"
    path "${sample_id}_filtered.fastq.gz"
    val outdir
    val gpu

    output:
    path "medaka_out/consensus.fasta"

    shell:
    """
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
    cpus = 4
    time = { 10.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path "medaka_out/consensus.fasta"
    path "${sample_id}_filtered.fastq.gz"
    val outdir

    output:
    path "${sample_id}_aligned.bam"
    path "${sample_id}_aligned.bam.bai"

    shell:
    """
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
    path "medaka_out/consensus.fasta"
    path "${sample_id}_aligned.bam"
    path "${sample_id}_aligned.bam.bai"
    val outdir

    output:
    path "${sample_id}_consensus_trimmed.fasta"
    path "${sample_id}_aligned.bam"
    path "${sample_id}_aligned.bam.bai"
    
    shell:
    """
    frags=(${fragments})

    for i in \${frags[@]}; do

        echo "Trimming \${i}..."
        python "${workflow.projectDir}/trim_contigs_by_depth.py" \
            -f "${outdir}/${sample_id}/\$i/${sample_id}_consensus.fasta" \
            -b "${outdir}/${sample_id}/\$i/${sample_id}_aligned.bam" \
            -o "${sample_id}_consensus_trimmed.fasta" \
            -d 1200
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
    cpus = 4
    time = { 10.minute * task.attempt}
    memory = { 16.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(fragments)
    path "${sample_id}_consensus_trimmed.fasta"
    path("${sample_id}_aligned.bam")
    path("${sample_id}_aligned.bam.bai")
    path reference
    path regions_bed
    val outdir
    val virus

    output:
    path "read_counts.txt"

    shell:
    """
    frags=(${fragments})
    for i in \${frags[@]}; do
        echo "Processing \${i}..."
#        if [ ${virus} == "SIV" ] && [ \$i == "Fragment_3" ]; then
            # Split reference into separate contigs
#            echo "Splitting Fragment 3"
#            seqtk subseq ${reference} ${regions_bed} > reference_contigs.fasta
#            python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}
#            for j in Fragment_3.1 Fragment_3.2; do
#                echo "Processing \${j}..."
#                echo "\$j" > frags.txt
#                seqtk subseq renamed_contigs.fasta frags.txt > ref.fasta
#
                # Count each 'haplotype'
#                python "${workflow.projectDir}/contig_read_counter.py" "${outdir}/${sample_id}/\$i/${sample_id}_aligned.bam" \
#                    --output-bam "${sample_id}_filtered.bam" \
#                    --count-output "read_counts_\${j}.txt"
#
#                # Get depth at each base along the consensus sequence
#                samtools sort "${sample_id}_filtered.bam" > "${sample_id}_filtered_sorted.bam"
#                samtools index "${sample_id}_filtered_sorted.bam"
#                samtools depth -a -@ 4 "${sample_id}_filtered_sorted.bam" > "${outdir}/${sample_id}/\$i/base_depth_\$j.txt"
#
#                cp "read_counts_\${j}.txt" "${outdir}/${sample_id}/\$i/"
#                # Add to make sure process completes correctly
#                cat "read_counts_\${j}.txt" >> read_counts.txt
#            done
#
#        else

            # Count each 'haplotype'
            python "${workflow.projectDir}/contig_read_counter.py" "${outdir}/${sample_id}/\$i/${sample_id}_aligned.bam" --output-bam "${sample_id}_filtered.bam"

            cp "read_counts.txt" "${outdir}/${sample_id}/\$i/"
#        fi
    done
    """

}

process firstConsensus {

    executor 'slurm'
    if(params.gpu) {
        conda "${workflow.projectDir}/rvhaplo.yaml"
        clusterOptions = '-A b1042 --gres=gpu:a100:1'
        queue = 'genomics-gpu'
    }
    else {
        conda "${workflow.projectDir}/rvhaplo-cpu.yaml"
        clusterOptions = '-A b1042'
        queue = 'genomics'
    }
    cpus = 4
    time = { 20.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2


    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
    path reference
    path regions_bed
    val virus
    path outdir
    val min_read_length
    val min_depth
    val gpu

    output:
    //path "${outdir}/${sample_id}/coverage.txt"
    path "${sample_id}*firstConsensus.fasta"

    shell:
    """

    # Get regions
    regions_array=(\$(awk '{ print \$4 }' $regions_bed))

    # Coverage of all regions
    minimap2 -ax lr:hq $reference ${sample_id}_concatenated.fastq.gz > full_genome.sam
    samtools view -e 'rlen>${min_read_length}' -bS -F 4 full_genome.sam > full_genome.bam
    samtools fastq full_genome.bam > ${outdir}/${sample_id}/${sample_id}_filtered.fastq
    samtools sort full_genome.bam > full_genome_sorted.bam
    samtools index full_genome_sorted.bam
    echo "Getting coverage"
    bedtools coverage -a $regions_bed -b full_genome_sorted.bam -nonamecheck > ${outdir}/${sample_id}/coverage.txt
    cat ${outdir}/${sample_id}/coverage.txt
    coverage_fractions=(\$(awk '{ print \$NF }' ${outdir}/${sample_id}/coverage.txt))
    coverage_depth=(\$(awk '{ print \$5 }' ${outdir}/${sample_id}/coverage.txt))

    # Split reference into separate contigs
    bedtools getfasta -fi ${reference} -bed ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}

    count_index=0
    for i in \${regions_array[@]}; do
        echo \${i}
        fraction=\${coverage_fractions[\$count_index]}
        fraction=\${fraction:0:4}
        fraction=\$(awk '{print \$1*100}' <<< \$fraction)
        depth=\${coverage_depth[\$count_index]}
        if [ "\$fraction" -lt "100" ] || [ "\$depth" -lt "$min_depth" ]; then
            let "count_index+=1"
            continue
        fi
        if [ ! -d ${outdir}/${sample_id}/\${i} ]; then
            mkdir ${outdir}/${sample_id}/\${i}
            echo "${sample_id},\$i" >> ${sample_id}Fragments.txt
        fi

        # Alignment and first consensus
        mini_align -i ${outdir}/${sample_id}/${sample_id}_filtered.fastq -r renamed_contigs.fasta -d lr:hq -m -t 4 -p firstAlign_\${i}
        medaka inference firstAlign_\${i}.bam ${sample_id}_\${i}.hdf --threads 2 --regions \${i} --model r1041_e82_400bps_hac_v4.3.0
        medaka sequence *_\${i}.hdf renamed_contigs.fasta ${sample_id}_\${i}_firstConsensus.fasta --threads 4 --no-fillgaps --regions \${i}
        cp ${sample_id}_\${i}_firstConsensus.fasta ${outdir}/${sample_id}/\${i}/${sample_id}_\${i}_firstConsensus.fasta

        # Error handling if consensus isn't created
        if [ ! -f ${outdir}/${sample_id}/\${i}/${sample_id}_\${i}_firstConsensus.fasta ]; then 
            echo "Could not build first consensus" >> $outdir/$sample_id/\${i}/${sample_id}_error.txt
            pwd >> ${outdir}/${sample_id}/\${i}/${sample_id}_error.txt
            echo "Could not build first consensus" > ${sample_id}_\${i}_firstConsensus.fasta
            continue
        fi
        let "count_index+=1"
    done

    """
}

process rvhaplo {
    
    executor 'slurm'
    if(params.gpu) {
        conda "${workflow.projectDir}/rvhaplo.yaml"
        clusterOptions = '-A b1042 --gres=gpu:a100:1'
        queue = 'genomics-gpu'
        cpus = 2
        time = { 60.minute * task.attempt }
        memory = { 40.GB * task.attempt }
    }
    else {
        conda "${workflow.projectDir}/rvhaplo-cpu.yaml"
        clusterOptions = '-A b1042'
        queue = 'genomics'
        cpus = 4
        time = { 120.minute * task.attempt }
        memory = { 70.GB * task.attempt }
    }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
    path "${sample_id}*firstConsensus.fasta"
    path outdir
    val subgraphs
    val abundance
    val smallest_snv
    path regions_bed
    val gpu

    shell:
    """

    # Get list of created fragment directories
    dir_list=(\$(ls ${outdir}/${sample_id}/ | grep Fragment))

    echo "\${dir_list[@]}"

    for i in "\${!dir_list[@]}"; do

        echo \${dir_list[\$i]}
        # Make sure first consensus exists
        if [ ! -f ${outdir}/${sample_id}/\${dir_list[\$i]}/*firstConsensus.fasta ]; then 
            echo "Could not build first consensus" > ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            pwd >> ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            echo "Could not build first consensus"
            continue
        fi

        # Test that fasta isn't just an error message
        line=\$(head -n 1 ${outdir}/${sample_id}/\${dir_list[\$i]}/*firstConsensus.fasta)
        if [ \$line == "Could not build first consensus" ]; then
            echo "Could not build first consensus" > ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            pwd >> ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            echo "Could not build first consensus"
            continue
        fi

        # Align to filter reads
        minimap2 -ax lr:hq ${outdir}/${sample_id}/\${dir_list[\$i]}/*firstConsensus.fasta ${outdir}/${sample_id}/${sample_id}_filtered.fastq > ${sample_id}_\${i}_IntermediateAlignment.sam
        # Remove unaligned reads and filter for length
        samtools view -bS -F 4 -O BAM -o ${sample_id}_\${i}_filtered.bam ${sample_id}_\${i}_IntermediateAlignment.sam
        samtools fastq ${sample_id}_\${i}_filtered.bam > ${sample_id}_\${i}_filtered.fastq
        cp ${sample_id}_\${i}_filtered.bam ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_filtered.bam
        cp ${sample_id}_\${i}_filtered.fastq ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_filtered.fastq

        # Rename contigs
        python ${workflow.projectDir}/rename_contigs.py -i ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}*firstConsensus.fasta -o ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta -b ${regions_bed}

        # Make sure fasta files have content
        firstcon_length=\$(wc -l < ${outdir}/${sample_id}/\${dir_list[\$i]}/*firstConsensus.fasta)
        filtered_length=\$(wc -l < ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_filtered.fastq)

        # Check that first consensus has content
        if [ \$firstcon_length -lt 2 ]; then
            echo "No reads aligned in first consensus" > ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_error.txt
        # Check that some filtered reads passed
        elif [ \$filtered_length -lt 2 ]; then
            echo "No reads aligned to region" > ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_error.txt
        else
            # Align to first consensus, make new consensus
            echo "Building consensus"
            mini_align -i ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_filtered.fastq -r ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta -d lr:hq -m -t 4 -p secondAlign
            medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
            medaka sequence ${sample_id}.hdf ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_realignment.fasta --threads 4
            minimap2 -ax lr:hq ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_realignment.fasta ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam
            rm ${outdir}/${sample_id}/\${dir_list[\$i]}/*.mmi

            # Run RVHaplo
            echo "Haplotype discovery"
            RVHaploPath=${workflow.projectDir}"/rvhaplo.sh"
            ("\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r ${outdir}/${sample_id}/\${dir_list[\$i]}/${sample_id}_realignment.fasta -o ${outdir}/${sample_id}/\${dir_list[\$i]} -p ${sample_id} -t 8 -e 0.1 -sg $subgraphs -a $abundance -ss $smallest_snv)

            echo "Finished "\${dir_list[\$i]}
        fi

    done

    echo "Finished haplotype discovery for all regions"

    """
}

process countBarcodes {

    executor 'slurm'
    if(params.gpu) {
        conda "${workflow.projectDir}/rvhaplo.yaml"
    }
    else {
        conda "${workflow.projectDir}/rvhaplo-cpu.yaml"
    }
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 2
    time = { 8.minute * task.attempt}
    memory = { 8.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
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
    perl ${workflow.projectDir}/bam_barcode_count.pl -b ${sample_id}_barcode_sorted.bam -s 113 -e 122 > ${outdir}/${sample_id}/barcode/barcodes.txt
    cp ${sample_id}_barcode_sorted.bam ${outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam
    cp ${sample_id}_barcode_sorted.bam.bai ${outdir}/${sample_id}/barcode/${sample_id}_barcode_sorted.bam.bai

    """
}

workflow {
    
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv()
        .set {barcodes_ch}

    concatenated = concatenateFastq(params.fastq_dir, barcodes_ch, params.outdir)
    quality = qualityFilter(barcodes_ch, concatenated, params.outdir)
    flye = flyeAssembly(barcodes_ch, quality, params.outdir)
    cleaned = cleanContigs(barcodes_ch, flye, params.reference, params.regions_bed, params.outdir)
    polished = polishAssembly(barcodes_ch, cleaned, quality, params.outdir, params.gpu)
    aligned = alignReads(barcodes_ch, polished, quality, params.outdir)
    trimmed = trimConsensus(barcodes_ch, polished, aligned, params.outdir)
    haplotypes(barcodes_ch, trimmed, params.reference, params.regions_bed, params.outdir, params.virus)

    //firstCon = firstConsensus(barcodes_ch, concatenated, params.reference, params.regions_bed, params.virus, params.outdir, params.min_read_length, params.min_depth, params.gpu)
    //haplotypes(barcodes_ch, concatenated, firstCon, params.outdir, params.subgraphs, params.abundance, params.smallest_snv, params.regions_bed, params.gpu)
    //if (params.split_barcode)
    //    countBarcodes(barcodes_ch, concatenated, params.outdir, params.gpu)
}
