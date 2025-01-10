# Haplotypes workflow

## Clone repository

```shell
git clone https://github.com/rlr-lab/haplotypes_workflow.git
cd haplotypes_workflow
```

## On Quest - activate Nextflow module

```shell
module purge all
# Version 24.04.4 is currently the most recent available on Quest
module load nextflow/24.04.4
nextflow run main.nf [options]
```

## Options

### General options

--barcodes &emsp; string, A 2-column csv file containing the fastq file prefixes (such as barcode01) and the sample name\
--fastq_dir &emsp; string, Path to a directory containing either a) a folder of fastq files for each barcode, or b) a fastq file for each barcode\
--regions_bed &emsp; string, Path to a .bed file with coordinates for the region of interest\
--reference &emsp; string, Path to a fasta file with a reference genome for reads to be aligned against (See Reference and Regions below)\
--outdir &emsp; string, Path for the output to be stored\
--virus &emsp; string, Either 'HIV' or 'SIV'\
-work-dir &emsp; string, Path to pipeline work directory (default: /projects/b1042/LorenzoRedondoLab/Seth/work)\
See more details [here](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters)

### Options passed to RV Haplo

--subgraphs &emsp; integer, Number of subgraphs to run MCL (default: 1)\
--abundance &emsp; float, A threshold for filtering low-abundance haplotypes. (default: 0.005)

## Reference and Regions

### HIV

Reference genome: HXB2Ref_FullGenome.fas\
Region .bed file: HXB2regions.bed

### SIV

Reference genome: SIVMac239FullGenome.fas\
Region .bed file: SIVregions.bed
