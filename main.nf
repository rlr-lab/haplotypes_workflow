params.barcodes = "barcodes.txt"
params.fastq_dir = ""
params.regions_bed = "${workflow.projectDir}/SIVregions.bed"
params.reference = "${workflow.projectDir}/SIVMac239FullGenome.fas"
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
    tuple val(barcode), val(sample_id), val(regions)
    val outdir
    errorStrategy 'retry'

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
    memory = { 2.GB * task.attempt}
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
    path "${sample_id}_firstConsensus.fasta"

    shell:
    """
    # Get regions
    regions_array=()
    IFS=" "
    if [ "${regions}" != "null" ]; then
        regions_array+=("${regions}")
    fi
    unset IFS

    # If no regions specified, do whole genome
    if [ -z "\${regions_array}" ]; then
        if [ "${virus}" == "HIV" ]; then
            regions_array+=("Fragment_1")
            regions_array+=("Fragment_2")
            regions_array+=("Fragment_3")
            regions_array+=("Fragment_4")
        elif [ "${virus}" == "SIV" ]; then
            regions_array+=("Fragment_1")
            regions_array+=("Fragment_2")
            regions_array+=("Fragment_3")
            regions_array+=("Fragment_4")
            regions_array+=("Fragment_5")
        fi
    fi

    # Split reference into separate contigs
    bedtools getfasta -fi ${reference} -bed ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}

    # Alignment and first consensus
    mini_align -i ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz -r renamed_contigs.fasta -m -t 4 -p firstAlign
    medaka inference firstAlign.bam ${sample_id}.hdf --threads 2 --regions ${regions} --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence *.hdf renamed_contigs.fasta ${sample_id}_firstConsensus.fasta --threads 4 --no-fillgaps --regions \${regions_array[@]}
    cp ${sample_id}_firstConsensus.fasta ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta

    # Error handling if consensus isn't created
    if [ ! -f ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta ]; then 
        echo "Could not build first consensus" >> $outdir/$sample_id/${sample_id}_error.txt
        pwd >> $outdir/$sample_id/${sample_id}_error.txt
        exit 1
    fi
    """
}

process haplotypes {
    
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
    tuple val(barcode), val(sample_id), val(regions)
    path "${sample_id}_concatenated.fastq.gz"
    path "${sample_id}_firstConsensus.fasta"
    path outdir
    val subgraphs
    val abundance
    val min_read_length
    val smallest_snv

    shell:
    """
    # Make sure first consensus exists
    if [ ! -f ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta ]; then 
        echo "Could not build first consensus" >> $outdir/$sample_id/${sample_id}_error.txt
        pwd >> $outdir/$sample_id/${sample_id}_error.txt
        exit 0
    fi

    # Align only to filter reads
    minimap2 -ax map-ont $outdir/$sample_id/${sample_id}_firstConsensus.fasta ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz > ${sample_id}_IntermediateAlignment.sam
    # Remove unaligned reads and filter for length
    samtools view -e 'rlen>${min_read_length}' -bS -F 4 -O BAM -o ${sample_id}_filtered.bam ${sample_id}_IntermediateAlignment.sam
    samtools fastq ${sample_id}_filtered.bam > ${sample_id}_filtered.fastq
    cp ${sample_id}_filtered.bam $outdir/$sample_id/
    cp ${sample_id}_filtered.fastq $outdir/$sample_id/

    # Rename contigs
    python ${workflow.projectDir}/rename_contigs.py -i ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta -o ${sample_id}_firstConsensus_renamed.fasta -c ${regions}

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
        mini_align -i ${sample_id}_filtered.fastq -r ${sample_id}_firstConsensus_renamed.fasta -m -t 4 -p secondAlign
        medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
        medaka sequence ${sample_id}.hdf ${sample_id}_firstConsensus_renamed.fasta $outdir/$sample_id/${sample_id}_realignment.fasta --threads 4
        minimap2 -ax map-ont $outdir/$sample_id/${sample_id}_realignment.fasta ${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam

        # Run RVHaplo
        RVHaploPath=${workflow.projectDir}"/rvhaplo.sh"
        . "\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r $outdir/$sample_id/${sample_id}_realignment.fasta \
		    -o $outdir/$sample_id -p ${sample_id} -t 8 -e 0.1 \
		    -sg $subgraphs -a $abundance -ss $smallest_snv;
    fi
    """
}

//process splitBarcode {
//
//    executor 'local'
//    conda "${workflow.projectDir}/rvhaplo.yaml"
//    errorStrategy 'retry'
//    maxRetries 2
//
//    input:
//    barcodes_ch
//    params.regions_bed
//
//    shell:
//    """
//    mkdir $outdir/$sample_id/barcodes
//    minimap2 -ax map-ont ${workflow.projectDir}/barcode_reference.fas ${sample_id}_concatenated.fastq.gz > ${sample_id}.sam
//    samtools view -F 4 -bS -h -O BAM -e 'rlen>233' ${sample_id}.sam > ${sample_id}.bam
//    samtools sort -o $outdir/$sample_id/barcodes/${sample_id}_sorted.bam ${sample_id}.bam
//    samtools index $outdir/$sample_id/barcodes/${sample_id}_sorted.bam
//    samtools depth -a -r "barcode" -@ 4 $outdir/$sample_id/barcodes/${sample_id}_barcode.bam > $outdir/$sample_id/barcodes/${sample_id}_barcode_depth.txt
//    bam trimBam $outdir/$sample_id/barcodes/${sample_id}_barcode.bam $outdir/$sample_id/barcodes/${sample_id}_ONLYbarcode.bam 112
//    samtools fasta -n -o $outdir/$sample_id/barcodes/${sample_id}_ONLYbarcode.fasta $outdir/$sample_id/barcodes/${sample_id}_ONLYbarcode.bam
//
//
//    samtools depth -a -r "barcode_0:1-34" -f ${sample_id}_barcodes.bam -@ 4 > $outdir/$sample_id/barcodes/${sample_id}_barcode_depth.txt
//
//
//    """
//}

workflow {
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv()
        .set {barcodes_ch}
    concatenated = concatenateFastq(params.fastq_dir, barcodes_ch, params.outdir)
    firstCon = firstConsensus(barcodes_ch, concatenated, params.reference, params.regions_bed, params.virus, params.outdir)
    haplotypes(barcodes_ch, concatenated, firstCon, params.outdir, params.subgraphs, params.abundance, params.min_read_length, params.smallest_snv)
}
