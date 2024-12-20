# Haplotypes workflow

## Clone repository

```shell
git clone https://github.com/rlr-lab/haplotypes_workflow.git
```

Or, pull to edit files before pushing to share

```shell
git pull https://github.com/rlr-lab/haplotypes_workflow.git
```

## Install conda environment

```shell
mamba env create -f rvhaplo.yaml
```

## On Quest - activate nextflow module and run

```shell
module load nextflow/24.04.4
nextflow run main.nf
```

## Reference and Regions

HXB2Ref_FullGenome.fas
: HIV reference genome
SIVMac239FullGenome.fas
: SIV reference genome

## Issues

Currently, most of the file paths are specific paths. These need to be updated to relative paths so the workflow can be run by another user.
