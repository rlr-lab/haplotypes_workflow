params.barcodes = "barcodes.txt"
params.fastq_dir = ""
params.regions_bed = "${workflow.projectDir}/ReferenceSequences/SIVregions_wBarcode.bed"
params.reference = "${workflow.projectDir}/ReferenceSequences/SIVMac239FullGenome_wBarcode.fas"
params.outdir = "${workflow.projectDir}/nf-results"
params.virus = "SIV"
params.min_read_length = 1200
params.min_depth = 1000
params.split_barcode = false
// RVHaplo parameters
params.subgraphs = 1
params.abundance = 0.001
params.smallest_snv = 5

process concatenateFastq {
    
    executor 'local'
    // Conda not needed, but initialized now so that the next process doesn't time out while creating environment
    //conda "${workflow.projectDir}/rvhaplo.yaml"
    conda "/home/lzh8485/.conda/envs/rvhaplo"

    input:
    path fastq_dir
    tuple val(barcode), val(sample_id)
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
        cat ${fastq_dir}/${barcode}/*.fastq.gz > ${sample_id}_concatenated.fastq.gz;
    elif [ -f ${fastq_dir}/${sample_id}.fastq.gz ]; then
        cat ${fastq_dir}/${sample_id}.fastq.gz > ${sample_id}_concatenated.fastq.gz;
    fi
    cp ${sample_id}_concatenated.fastq.gz ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz

    """

}

process firstConsensus {

    executor 'slurm'
    //conda "${workflow.projectDir}/rvhaplo.yaml"
    conda "/home/lzh8485/.conda/envs/rvhaplo"
    clusterOptions = '-A b1042'
    queue = 'genomics'
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

    output:
    //path "${outdir}/${sample_id}/coverage.txt"
    path "${sample_id}*firstConsensus.fasta"

    shell:
    """

    # Get regions
    regions_array=(\$(awk '{ print \$4 }' $regions_bed))

    # Coverage of all regions
    minimap2 -ax lr:hq $reference ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz > full_genome.sam
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
            continue
        fi
        let "count_index+=1"
    done

    """
}

process haplotypes {
    
    executor 'slurm'
    //conda "${workflow.projectDir}/rvhaplo.yaml"
    conda "/home/lzh8485/.conda/envs/rvhaplo"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 8
    time = { 120.minute * task.attempt }
    memory = { 60.GB + (60.GB * task.attempt) }
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

    shell:
    """

    # Get list of created fragment directories
    dir_list=(\$(echo ${outdir}/${sample_id}/Fragment*))

    echo "\${dir_list[@]}"

    for i in "\${!dir_list[@]}"; do

        echo \${dir_list[\$i]}
        # Make sure first consensus exists
        if [ \${dir_list[\$i]} == "${outdir}/${sample_id}/barcode" ]; then
            echo "barcode dir"
            continue
        elif [ ! -f \${dir_list[\$i]}/*firstConsensus.fasta ]; then 
            echo "Could not build first consensus" > \${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            pwd >> \${dir_list[\$i]}/${sample_id}_\${i}_error.txt
            echo "Could not build first consensus"
            continue
        fi

        # Align to filter reads
        minimap2 -ax lr:hq \${dir_list[\$i]}/*firstConsensus.fasta ${outdir}/${sample_id}/${sample_id}_filtered.fastq > ${sample_id}_\${i}_IntermediateAlignment.sam
        # Remove unaligned reads and filter for length
        samtools view -bS -F 4 -O BAM -o ${sample_id}_\${i}_filtered.bam ${sample_id}_\${i}_IntermediateAlignment.sam
        samtools fastq ${sample_id}_\${i}_filtered.bam > ${sample_id}_\${i}_filtered.fastq
        cp ${sample_id}_\${i}_filtered.bam \${dir_list[\$i]}/${sample_id}_filtered.bam
        cp ${sample_id}_\${i}_filtered.fastq \${dir_list[\$i]}/${sample_id}_filtered.fastq

        # Rename contigs
        python ${workflow.projectDir}/rename_contigs.py -i \${dir_list[\$i]}/${sample_id}*firstConsensus.fasta -o \${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta -b ${regions_bed}

        # Make sure fasta files have content
        firstcon_length=\$(wc -l < \${dir_list[\$i]}/*firstConsensus.fasta)
        filtered_length=\$(wc -l < \${dir_list[\$i]}/${sample_id}_filtered.fastq)

        # Check that first consensus has content
        if [ \$firstcon_length -lt 2 ]; then
            echo "No reads aligned in first consensus" > \${dir_list[\$i]}/${sample_id}_error.txt
        # Check that some filtered reads passed
        elif [ \$filtered_length -lt 2 ]; then
            echo "No reads aligned to region" > \${dir_list[\$i]}/${sample_id}_error.txt
        else
            # Align to first consensus, make new consensus
            echo "Building consensus"
            mini_align -i \${dir_list[\$i]}/${sample_id}_filtered.fastq -r \${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta -d lr:hq -m -t 4 -p secondAlign
            medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
            medaka sequence ${sample_id}.hdf \${dir_list[\$i]}/${sample_id}_firstConsensus_renamed.fasta \${dir_list[\$i]}/${sample_id}_realignment.fasta --threads 4
            minimap2 -ax lr:hq \${dir_list[\$i]}/${sample_id}_realignment.fasta \${dir_list[\$i]}/${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam

            # Run RVHaplo
            echo "Haplotype discovery"
            RVHaploPath=${workflow.projectDir}"/rvhaplo.sh"
            ("\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r \${dir_list[\$i]}/${sample_id}_realignment.fasta -o \${dir_list[\$i]} -p ${sample_id} -t 8 -e 0.1 -sg $subgraphs -a $abundance -ss $smallest_snv)

            echo "Finished "\${dir_list[\$i]}
        fi

    done

    echo "Finished haplotype discovery for all regions"

    """
}

process countBarcodes {

    executor 'slurm'
    //conda "${workflow.projectDir}/rvhaplo.yaml"
    conda "/home/lzh8485/.conda/envs/rvhaplo"
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

    shell:
    """

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
    
    def bedfile = new File(params.regions_bed)
    def bed_rows = bedfile.readLines()*.split('	')
    Set regions = bed_rows*.getAt(3)
    Channel
        .fromList(regions)
        .set {fragment_ch}


    def barcode_file = new File(params.barcodes)
    def barcode_rows = barcode_file.readLines()*.split(',')
    Set barcode_name = barcode_rows*.getAt(0)
    Channel
        .fromList(barcode_name)
        .set {barcode_ch}

    Set sample_names = barcode_rows*.getAt([0,1])
    Channel
        .fromList(sample_names)
        .set {sample_ch}

    barcode_ch.join(sample_ch).set {rename_ch}
    barcode_ch.join(sample_ch).combine(fragment_ch).transpose(remainder: true).set {final_fragment_ch}
    

    concatenated = concatenateFastq(params.fastq_dir, rename_ch, params.outdir)
    firstCon = firstConsensus(rename_ch, concatenated, params.reference, params.regions_bed, params.virus, params.outdir, params.min_read_length, params.min_depth)
    haplotypes(rename_ch, concatenated, firstCon, params.outdir, params.subgraphs, params.abundance, params.smallest_snv, params.regions_bed)
    if (params.split_barcode)
        countBarcodes(rename_ch, concatenated, params.outdir)
}
