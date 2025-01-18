params.barcodes = "/projects/b1042/LorenzoRedondoLab/Seth/F3_bams/barcodes.txt"
params.fastq_dir = "/projects/b1042/LorenzoRedondoLab/Seth/F3_bams"
params.regions_bed = "/home/lzh8485/haplotypes_workflow/SIVregions.bed"
params.reference = "/home/lzh8485/haplotypes_workflow/SIVMac239FullGenome.fas"
params.outdir = "/projects/b1042/LorenzoRedondoLab/Seth/F3_bams/nf-results"
params.virus = "SIV"
// Match RV Haplo defaults but allow changing
params.subgraphs = 1
params.abundance = 0.01

process concatenateFastq {
    
    executor 'local'
    // Conda not needed, but initialized now so that the next process doesn't time out while creating environment
    conda 'rvhaplo.yaml'

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
    elif [ d ${outdir}/${sample_id} ]; then
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
    conda 'rvhaplo.yaml'
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 4
    time = { 5.minute * task.attempt}
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
    regions_array+=("${regions}")
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
            regions_array+=("Fragment_3.1")
            regions_array+=("Fragment_3.2")
            regions_array+=("Fragment_4")
            regions_array+=("Fragment_5")
        fi
    # Split SIV Fragment_3 into 3 parts
    else
        for i in \${regions_array}; do
            if [ "\$i" == "Fragment_3" ] && [ "${virus}" == "SIV" ]; then
                todelete="Fragment_3"
                regions_array=( "\$regions_array[@]/\$todelete}" )
                regions_array+=("Fragment_3.1")
                regions_array+=("Fragment_3.2")
                regions_array+=("barcode")
            fi
        done
    fi

    # Split reference into separate contigs
    bedtools getfasta -fi ${reference} -bed ${regions_bed} > reference_contigs.fasta
    python ${workflow.projectDir}/rename_contigs.py -i reference_contigs.fasta -o renamed_contigs.fasta -b ${regions_bed}

    # Alignment and first consensus
    mini_align -i ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz -r renamed_contigs.fasta -m -t 4 -p firstAlign
    medaka inference firstAlign.bam ${sample_id}.hdf --threads 2 --regions ${regions} --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence *.hdf renamed_contigs.fasta ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta --threads 4 --no-fillgaps --regions ${regions}
    cp ${outdir}/${sample_id}/${sample_id}_firstConsensus.fasta ${sample_id}_firstConsensus.fasta
    """
}

process haplotypes {
    
    executor 'slurm'
    conda 'rvhaplo.yaml'
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 8
    time = { 120.minute * task.attempt}
    memory = { 80.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id), val(regions)
    path "${sample_id}_concatenated.fastq.gz"
    path "${sample_id}_firstConsensus.fasta"
    path outdir
    val subgraphs
    val abundance

    shell:
    """
    # Remove sequences <1800b
    bioawk -c fastx '(length(\$seq)>1800) {print ">" \$name ORS \$seq}' ${sample_id}_concatenated.fastq.gz > ${sample_id}_filtered.fastq

    # Align to first consensus, make new consensus
    mini_align -i ${sample_id}_filtered.fastq -r ${sample_id}_firstConsensus.fasta -m -t 4 -p secondAlign
    medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence ${sample_id}.hdf ${sample_id}_firstConsensus.fasta $outdir/$sample_id/${sample_id}_realignment.fasta --threads 4
    minimap2 -ax map-ont $outdir/$sample_id/${sample_id}_realignment.fasta ${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam

    # Run RVHaplo
    RVHaploPath=${workflow.projectDir}"/rvhaplo.sh"
    . "\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r $outdir/$sample_id/${sample_id}_realignment.fasta -o $outdir/$sample_id -p ${sample_id} -t 8 -sg $subgraphs -a $abundance;
    """
}

workflow {
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv()
        .set {barcodes_ch}
    concatenated = concatenateFastq(params.fastq_dir, barcodes_ch, params.outdir)
    firstCon = firstConsensus(barcodes_ch, concatenated, params.reference, params.regions_bed, params.virus, params.outdir)
    haplotypes(barcodes_ch, concatenated, firstCon, params.outdir, params.subgraphs, params.abundance)
}
