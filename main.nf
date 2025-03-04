params.barcodes = "barcodes.txt"
params.fastq_dir = ""
params.regions_bed = "${workflow.projectDir}/ReferenceSequences/SIVregions_wBarcode.bed"
params.reference = "${workflow.projectDir}/ReferenceSequences/SIVMac239FullGenome_wBarcode.fas"
params.outdir = "${workflow.projectDir}/nf-results"
params.virus = "SIV"
params.min_read_length = 1200
params.split_barcode = false
// RVHaplo parameters
params.subgraphs = 1
params.abundance = 0.001
params.smallest_snv = 5

process concatenateFastq {
    
    executor 'local'
    // Conda not needed, but initialized now so that the next process doesn't time out while creating environment
    conda "${workflow.projectDir}/rvhaplo.yaml"

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
    conda "${workflow.projectDir}/rvhaplo.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 4
    time = { 8.minute * task.attempt}
    memory = { 4.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2


    input:
    tuple val(barcode), val(sample_id), val(regions)
    path "${sample_id}_concatenated.fastq.gz"
    path reference
    path regions_bed
    val virus
    path outdir

    output:
    //path "${outdir}/${sample_id}/coverage.txt"

    shell:
    """
    # Get regions
    regions_array=(\$(awk '{ print \$4 }' $regions_bed))

    # Coverage of all regions
    minimap2 -ax lr:hq $reference ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz > full_genome.sam
    samtools view -bS -F 4 full_genome.sam > full_genome.bam
    samtools sort full_genome.bam > full_genome_sorted.bam
    samtools index full_genome_sorted.bam
    bedtools coverage -a $regions_bed -b full_genome_sorted.bam -nonamecheck > ${outdir}/${sample_id}/coverage.txt
    coverage_fractions=(\$(awk '{ print \$NF }' ${outdir}/${sample_id}/coverage.txt))

    # Split reference into separate contigs
    bedtools getfasta -fi ${reference} -bed ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}

    count_index=0
    for i in \${regions_array[@]}; do
        fraction=\${coverage_fractions[\$count_index]}
        fraction=\${fraction:0:4}
        fraction=\$(awk '{print \$1*100}' <<< \$fraction)
        if [ "\$fraction" -lt "100" ]; then
            !((count_index++))
            continue
        fi
        if [ ! -d ${outdir}/${sample_id}/\${i} ]; then
            mkdir ${outdir}/${sample_id}/\${i}
            echo "${sample_id},\$i" >> ${sample_id}Fragments.txt
        fi

        # Alignment and first consensus
        mini_align -i ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz -r renamed_contigs.fasta -d lr:hq -m -t 4 -p firstAlign_\${i}
        medaka inference firstAlign_\${i}.bam ${sample_id}_\${i}.hdf --threads 2 --regions \${i} --model r1041_e82_400bps_hac_v4.3.0
        medaka sequence *_\${i}.hdf renamed_contigs.fasta ${sample_id}_\${i}_firstConsensus.fasta --threads 4 --no-fillgaps --regions \${i}
        cp ${sample_id}_\${i}_firstConsensus.fasta ${outdir}/${sample_id}/\${i}/${sample_id}_\${i}_firstConsensus.fasta

        # Error handling if consensus isn't created
        if [ ! -f ${outdir}/${sample_id}/\${i}/${sample_id}_\${i}_firstConsensus.fasta ]; then 
            echo "Could not build first consensus" >> $outdir/$sample_id/\${i}/${sample_id}_error.txt
            pwd >> ${outdir}/${sample_id}/\${i}/${sample_id}_error.txt
            continue
        fi
        !((count_index++))
    done
    """
}

process haplotypes {
    
    tag { fragment }

    executor 'slurm'
    conda "${workflow.projectDir}/rvhaplo.yaml"
    clusterOptions = '-A b1042'
    queue = 'genomics-himem'
    cpus = 8
    time = { 180.minute + (300.minute * (task.attempt - 1)) }
    memory = { 180.GB + (60.GB * task.attempt) }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
    val fragment
    path outdir
    val subgraphs
    val abundance
    val min_read_length
    val smallest_snv

    shell:
    """

    # Get list of created fragment directories
    dir_list=(\$(echo ${outdir}/${sample_id}/*/))

    count_index=0
    for i in \${dir_list[@]}; do

        # Make sure first consensus exists
        if [ ! -f \${i}${sample_id}_*_firstConsensus.fasta ]; then 
            echo "Could not build first consensus" >> $outdir/$sample_id/\${i}${sample_id}_\${count_index}_error.txt
            pwd >> \${i}${sample_id}_\${count_index}_error.txt
            !((count_index++))
            continue
        fi

        # Align only to filter reads
        minimap2 -ax lr:hq \${i}${sample_id}_*_firstConsensus.fasta ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz > ${sample_id}_\${count_index}_IntermediateAlignment.sam
        # Remove unaligned reads and filter for length
        samtools view -e 'rlen>${min_read_length}' -bS -F 4 -O BAM -o ${sample_id}_\${count_index}_filtered.bam ${sample_id}_\${count_index}_IntermediateAlignment.sam
        samtools fastq ${sample_id}_\${count_index}_filtered.bam > ${sample_id}_\${count_index}_filtered.fastq
        cp ${sample_id}_\${count_index}_filtered.bam \${i}${sample_id}_filtered.bam
        cp ${sample_id}_\${count_index}_filtered.fastq \${i}${sample_id}_filtered.fastq

        # Rename contigs
        python ${workflow.projectDir}/rename_contigs.py -i ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta -o ${sample_id}_firstConsensus_renamed.fasta -b ${regions_bed}

        # Make sure fasta files have content
        firstcon_length=\$(wc -l < ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta)
        filtered_length=\$(wc -l < ${sample_id}_filtered.fastq)

        # Check that first consensus has content
        if [ \$firstcon_length -lt 2 ]; then
            echo "No reads aligned in first consensus" >> $outdir/$sample_id/${sample_id}_error.txt
        # Check that some filtered reads passed
        elif [ \$filtered_length -lt 2 ]; then
            echo "QC Fail - Possibly no long reads (length > ${min_read_length}) present" >> $outdir/$sample_id/${sample_id}_error.txt
        else
            # Align to first consensus, make new consensus
            mini_align -i ${sample_id}_filtered.fastq -r ${sample_id}_firstConsensus_renamed.fasta -d lr:hq -m -t 4 -p secondAlign
            medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
            medaka sequence ${sample_id}.hdf ${sample_id}_firstConsensus_renamed.fasta $outdir/$sample_id/${sample_id}_realignment.fasta --threads 4
            minimap2 -ax lr:hq $outdir/$sample_id/${sample_id}_realignment.fasta ${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam

            # Run RVHaplo
            RVHaploPath=${workflow.projectDir}"/rvhaplo.sh"
            . "\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r $outdir/$sample_id/${sample_id}_realignment.fasta \
		        -o $outdir/$sample_id -p ${sample_id} -t 8 -e 0.1 \
		        -sg $subgraphs -a $abundance -ss $smallest_snv;
        fi
    
    """
}

process countBarcodes {

    executor 'slurm'
    conda "${workflow.projectDir}/rvhaplo.yaml"
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

    minimap2 -ax lr:hq ${workflow.projectDir}/barcode_reference.fas ${sample_id}_concatenated.fastq.gz > ${sample_id}_barcode.sam
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
    firstCon = firstConsensus(rename_ch, concatenated, params.reference, params.regions_bed, params.virus, params.outdir)
    //haplotypes(barcodes_ch, concatenated, fragment_ch, params.outdir, params.subgraphs, params.abundance, params.min_read_length, params.smallest_snv)
    if (params.split_barcode)
        countBarcodes(barcodes_ch, concatenated, params.outdir)
}
