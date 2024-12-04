params.barcodes = "/projects/b1042/LorenzoRedondoLab/11272024_EPG_PBMCs_Tissue/NanoporeBarcodesFormat.csv"
params.fastq_dir = "/projects/b1042/LorenzoRedondoLab/11272024_EPG_PBMCs_Tissue/11272024_EPG_PBMCs_Tissue/20241127_1916_X4_FAZ62876_9d2b31c7/fastq_pass"
params.regions_bed = "/home/lzh8485/haplotypes_workflow/HXB2regions.bed"
params.reference = "/home/lzh8485/haplotypes_workflow/HXB2Ref_FullGenome.fas"
params.outdir = "/projects/b1042/LorenzoRedondoLab/11272024_EPG_PBMCs_Tissue/nf-results"

process concatenateFastq {
    
    executor 'local'

    input:
    path fastq_dir
    tuple val(barcode), val(sample_id)
    val outdir
    errorStrategy 'retry'

    output:
    path "${sample_id}_concatenated.fastq.gz"

    shell:
    """
    if [ ! -d $outdir]; then
        mkdir $outdir
    fi
    if [ ! -d $outdir/$sample_id ]; then
        mkdir $outdir/$sample_id
    elif [ d $outdir/$sample_id ]; then
        rm -rf $outdir/$sample_id
        mkdir $outdir/$sample_id
    fi
    if [ -d $fastq_dir/$barcode ]; then
        cat $fastq_dir/$barcode/*.fastq.gz > ${sample_id}_concatenated.fastq.gz
    elif [ -f $fastq_dir/${sample_id}.fastq.gz ]; then
        cat $fastq_dir/${sample_id}.fastq.gz > ${sample_id}_concatenated.fastq.gz
    fi
    cp ${sample_id}_concatenated.fastq.gz ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz
    """

}

process getRegions {

    executor 'local'
    
    input:
    tuple val(barcode), val(sample_id)
    path regions_bed
    path fastq_dir
    errorStrategy 'retry'

    output:
    env 'filter_pos'

    shell:
    """
    if [ "$fastq_dir" == "F1_bams" ]; then
        filter_pos=\$(awk -v pat="Fragment_1" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    elif [ "$fastq_dir" == "F2_bams" ]; then
        filter_pos=\$(awk -v pat="Fragment_2" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    elif [ "$fastq_dir" == "F3_bams" ]; then
        filter_pos=\$(awk -v pat="Fragment_3" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    elif [ "$fastq_dir" == "F4_bams" ]; then
        filter_pos=\$(awk -v pat="Fragment_4" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    elif [ "$fastq_dir" == "F5_bams" ]; then
        filter_pos=\$(awk -v pat="Fragment_5" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    else
        filter_pos=\$(awk -v pat="Fragment_1" '\$4 ~ pat {print \$1 ":" \$2 "-" \$3}' $regions_bed)
    fi
    """

}

process firstConsensus {

    executor 'slurm'
    conda '/home/lzh8485/.conda/envs/rvhaplo'
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 4
    time = { 5.minute * task.attempt}
    memory = { 2.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2


    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
    val filter_pos
    path reference
    path outdir

    output:
    path "${sample_id}_firstConsensus.fasta"

    shell:
    """
    mini_align -i ${outdir}/${sample_id}/${sample_id}_concatenated.fastq.gz -r $reference -m -t 4 -p firstAlign
    medaka inference firstAlign.bam $sample_id\.hdf --threads 2 --regions $filter_pos --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence ${sample_id}.hdf $reference $outdir/$sample_id/${sample_id}_firstConsensus.fasta --threads 4 --no-fillgaps --regions $filter_pos
    cp $outdir/$sample_id/${sample_id}_firstConsensus.fasta ${sample_id}_firstConsensus.fasta
    """
}

process haplotypes {
    
    executor 'slurm'
    conda '/home/lzh8485/.conda/envs/rvhaplo'
    clusterOptions = '-A b1042'
    queue = 'genomics'
    cpus = 8
    time = { 30.minute * task.attempt}
    memory = { 80.GB * task.attempt}
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple val(barcode), val(sample_id)
    path "${sample_id}_concatenated.fastq.gz"
    path "${sample_id}_firstConsensus.fasta"
    path outdir

    shell:
    """
    bioawk -c fastx '(length(\$seq)>1800) {print ">" \$name ORS \$seq}' ${sample_id}_concatenated.fastq.gz > ${sample_id}_filtered.fastq
    mini_align -i ${sample_id}_filtered.fastq -r ${sample_id}_firstConsensus.fasta -m -t 4 -p secondAlign
    medaka inference secondAlign.bam ${sample_id}.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence ${sample_id}.hdf ${sample_id}_firstConsensus.fasta $outdir/$sample_id/${sample_id}_realignment.fasta --threads 4
    minimap2 -ax map-ont $outdir/$sample_id/${sample_id}_realignment.fasta ${sample_id}_filtered.fastq > ${sample_id}_RVHaploinput.sam
    RVHaploPath="/home/lzh8485/haplotypes_workflow/rvhaplo.sh";
    . "\${RVHaploPath}" -i ${sample_id}_RVHaploinput.sam -r $outdir/$sample_id/${sample_id}_realignment.fasta -o $outdir/$sample_id -p ${sample_id} -t 8;
    """
}

workflow {
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv()
        .set {barcodes_ch}
    concatenated = concatenateFastq(params.fastq_dir, barcodes_ch, params.outdir)
    regions = getRegions(barcodes_ch, params.regions_bed, params.fastq_dir)
    firstCon = firstConsensus(barcodes_ch, concatenated, regions, params.reference, params.outdir)
    haplotypes(barcodes_ch, concatenated, firstCon, params.outdir)
}