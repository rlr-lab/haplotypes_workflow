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

|Option|Description|
|:-------|:-----------|
|--barcode|*string*, A 3-column csv file containing the fastq file prefixes (such as barcode01), the sample name, and the fragment(s) sequenced. See example below|
|--fastq_dir|*string*, Path to a directory containing either a) a folder of fastq files for each barcode, or b) a fastq file for each barcode|
|--regions_bed|*string*, Path to a .bed file with coordinates for the region of interest|
|--reference|*string*, Path to a fasta file with a reference genome for reads to be aligned against (See Reference and Regions below)|
|--outdir|*string*, Path for the output to be stored|
|--virus|*string*, Either 'HIV' or 'SIV'|
|-work-dir|*string*, Path to pipeline work directory (default: /projects/b1042/LorenzoRedondoLab/Seth/work)|

See more details [here](https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters)

### Options passed to RV Haplo

|Option|Description|
|:-----|:----------|
|--subgraphs|*integer*, Number of subgraphs to run MCL (default: 1)|
|--abundance|*float*, A threshold for filtering low-abundance haplotypes. (default: 0.005)|

## Reference and Regions

### HIV

Reference genome: HXB2Ref_FullGenome.fas\
Region .bed file: HXB2regions.bed

### SIV

Reference genome: SIVMac239FullGenome.fas\
Region .bed file: SIVregions.bed

### Barcode .csv Example

```text
barcode01,SampleA,Fragment_1
barcode02,SampleA,Fragment_2
barcode03,SampleB,Fragment_4 Fragment_5
```
