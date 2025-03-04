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

An alternative option is to use a conda environment instead of a pre-installed module. This allows the program to run on other HPC systems besides Northwestern's Quest HPC.

```shell
mamba create -n nextflow nextflow
mamba activate nextflow
```

If you are submitting a job to a queue to run the pipeline, add the following to the beginning of your workflow to activate your new conda environment.

```shell
module load mamba
source /home/[USERID]/.bashrc
conda activate /home/[USERID]/.conda/envs/nextflow
```

## Options

### General options

|Option|Description|
|:-------|:-----------|
|--barcodes|*string*, A 2-column csv file containing the fastq file prefixes (such as barcode01) and the sample name. See example below|
|--fastq_dir|*string*, Path to a directory containing either (a) a folder of fastq files for each barcode, or (b) a fastq file for each barcode|
|--regions_bed|*string*, Path to a .bed file with coordinates for the region of interest (default: SIVregions.bed)|
|--reference|*string*, Path to a fasta file with a reference genome for reads to be aligned against. See Reference and Regions below (default: SIVMac239FullGenome.fas)|
|--outdir|*string*, Path for the output to be stored (default: nf-results/)|
|--virus|*string*, Either 'HIV' or 'SIV' (default: SIV)|
|--min_read_length|*integer*, Minimum read length allowable for QC filtering (default: 1200)|
|--split_barcode|*boolean*, Should the barcode analysis for SIV Fragment 3 be performed (default: false)|
|-work-dir|*string*, Path to pipeline work directory (default: /projects/b1042/LorenzoRedondoLab/Seth/work)|

See more details [here](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters)

### Options passed to RV Haplo

|Option|Description|
|:-----|:----------|
|--subgraphs|*integer*, Number of subgraphs to run MCL (default: 1)|
|--abundance|*float*, A threshold for filtering low-abundance haplotypes. (default: 0.001)|
|--smallest_snv|*integer*, Minimum # of SNV sites for haplotype construction. (default: 5)|

See the RVHaplo documentation [here](https://github.com/dhcai21/RVHaplo) and the journal article describing the software [here](https://doi.org/10.1093/bioinformatics/btac089).

## Reference and Regions

Found in the `ReferenceSequences/` folder

### HIV

Reference genome: HXB2Ref_FullGenome.fas\
Region .bed file: HXB2regions.bed

### SIV

Reference genome: SIVMac239FullGenome_wBarcode.fas\
Region .bed file: SIVregions_wBarcode.bed

### Barcode .csv Example

```text
barcode01,SampleA
barcode02,SampleB
barcode03,SampleC
```
