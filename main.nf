params.barcodes = "/projects/p32511/F2_bams/barcodes.txt"
params.fastq_dir = "/projects/p32511/F2_bams"
params.regions_bed = "/home/lzh8485/haplotypes_workflow/SIVregions.bed"
params.reference = "/home/lzh8485/haplotypes_workflow/SIVMac239FullGenome.fas"
params.outdir = "/projects/p32511/F2_bams/nf-results"

process concatenateFastq {
    
    executor 'local'

    input:
    path fastq_dir
    tuple val(barcode), val(sample_id)

    output:
    path "concatenated.fastq.gz"

    shell:
    """
    if [ -d "$fastq_dir/$barcode" ]; then
        cat $fastq_dir/$barcode/*.fastq.gz > concatenated.fastq.gz
    elif [ -f "$fastq_dir/${sample_id}.fastq.gz" ]; then
        cat $fastq_dir/${sample_id}.fastq.gz > concatenated.fastq.gz
    fi
    """

}

process getRegions {

    executor 'local'
    
    input:
    tuple val(barcode), val(sample_id)
    path regions_bed
    path fastq_dir

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
    path 'concatenated.fastq.gz'
    val filter_pos
    path reference
    path outdir

    output:
    path 'firstConsensus.fasta'

    shell:
    """
    if [ ! -d "$outdir"]; then
        mkdir $outdir
    fi
    if [ ! -d "$outdir/$sample_id" ]; then
        mkdir $outdir/$sample_id
    fi
    mini_align -i concatenated.fastq.gz -r $reference -m -t 4 -p firstAlign
    medaka inference firstAlign.bam $sample_id\.hdf --threads 2 --regions $filter_pos --model r1041_e82_400bps_hac_v4.3.0
    medaka sequence ${sample_id}.hdf $reference $outdir/$sample_id/firstConsensus.fasta --threads 4
    cp $outdir/$sample_id/firstConsensus.fasta firstConsensus.fasta
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
    path 'concatenated.fastq.gz'
    path 'firstConsensus.fasta'
    path outdir

    shell:
    """
    if [ ! -f $outdir/$sample_id/rvhaplo_haplotypes.fasta ]; then    
        minimap2 -ax map-ont firstConsensus.fasta concatenated.fastq.gz > realignment.sam;
        RVHaploPath="/home/lzh8485/haplotypes_workflow/rvhaplo.sh";
        . "\${RVHaploPath}" -i realignment.sam -r firstConsensus.fasta -o $outdir/$sample_id -t 8 -a 0.01 -sg 32;
    fi
    """
}

workflow {
    Channel
        .fromPath(params.barcodes, checkIfExists: true, type: 'file')
        .splitCsv()
        .set {barcodes_ch}
    concatenated = concatenateFastq(params.fastq_dir, barcodes_ch)
    regions = getRegions(barcodes_ch, params.regions_bed, params.fastq_dir)
    firstCon = firstConsensus(barcodes_ch, concatenated, regions, params.reference, params.outdir)
    haplotypes(barcodes_ch, concatenated, firstCon, params.outdir)
}