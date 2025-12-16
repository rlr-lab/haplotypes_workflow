# Long-Read Viral Consensus Pipeline (Nextflow DSL2)

This repository contains a Nextflow DSL2 pipeline for generating viral consensus sequences from long-read sequencing data (e.g. Oxford Nanopore). It supports both de novo assembly and reference-guided assembly, with optional polishing and read realignment.

The pipeline is designed for HPC environments and is compatible with Singularity/Apptainer containers.

---

## Overview

For each sample, the pipeline performs:

1. FASTQ concatenation
2. Read quality and length filtering
3. Either:
   - De novo assembly with Flye → Medaka polishing → read realignment, or
   - Reference-guided assembly with minimap2 + iVar consensus
4. Output of final BAM files and consensus FASTA sequences

The workflow automatically selects de novo mode if no reference is provided.

---

## Workflow Summary

### Common steps

- Concatenate FASTQ files per sample
- Quality and length filtering using fastplong

### De novo mode (`--denovo true` or no reference)

1. Assembly with Flye
2. Polishing with Medaka
3. Read alignment to polished consensus

### Reference-guided mode (`--reference reference.fasta`)

1. Initial alignment to reference
2. First consensus with iVar
3. Realignment to first consensus
4. Final consensus with iVar

---

## Input Requirements

### 1. Barcode/sample CSV

A two-column CSV file (no header):

```{text}
barcode01,sampleA
barcode02,sampleB
```

- Column 1: barcode or FASTQ directory name
- Column 2: sample ID (used for output naming)

Provide via:

```{shell}
--barcodes samples.csv
```

### 2. FASTQ layout

Either of the following structures is supported:

#### Option A: Per-barcode directories

```{shell}
fastq/
├── barcode01/
│   ├── read1.fastq.gz
│   ├── read2.fastq.gz
```

#### Option B: Single FASTQ per sample

```{shell}
fastq/
├── sampleA.fastq.gz
```

Provide via:

```{shell}
--fastq_dir /path/to/fastq
```

---

## Output Structure

All results are written to:

```{shell}
--outdir nf-results/
```

Per-sample layout:

```{shell}
nf-results/
└── sampleA/
    ├── fastplong/
    │   ├── sampleA_filtered.fastq.gz
    │   ├── sampleA_fastp_report.html
    │   └── sampleA_fastp_report.json
    ├── flye_out/             # (de novo only)
    ├── sampleA_consensus.fasta
    ├── sampleA_aligned.bam
    └── sampleA_aligned.bam.bai
```

---

## Parameters

| Parameter           | Description                        | Default      |
| ------------------- | ---------------------------------- | ------------ |
| `--barcodes`        | CSV mapping barcodes to sample IDs | **required** |
| `--fastq_dir`       | Directory containing FASTQ files   | **required** |
| `--outdir`          | Output directory                   | `nf-results` |
| `--reference`       | Reference FASTA (optional)         | none         |
| `--denovo`          | Force de novo assembly             | auto         |
| `--min_read_length` | Minimum read length                | 1200         |
| `--max_read_length` | Maximum read length                | 20000        |

If `--reference` is not provided, the pipeline automatically runs in de novo mode.

---

## This pipeline assumes execution inside a container

There is an included container that includes:

- fastplong
- flye
- medaka
- minimap2
- samtools
- ivar
- seqtk
- Python 3

A prebuilt OCI-compatible container (Docker → Singularity) is recommended.

---

## Running the Pipeline

Example:

```{shell}
nextflow run main.nf \
  --barcodes samples.csv \
  --fastq_dir fastq \
  --outdir nf-results \
  -work /tmp/work
```

---

## Notes and Assumptions

- The pipeline assumes one viral genome per sample
- Flye failures fall back to using the provided reference (if available)
- Read filtering is tuned for long-read viral sequencing and may require adjustment for other use cases
