# SIV Tiled Amplicon Assembly Pipeline

This Nextflow pipeline processes SIV tiled amplicon sequencing data from raw FASTQ files to high-quality consensus sequences. It includes quality filtering, de novo assembly, polishing, trimming, haplotype reconstruction, and optional barcode counting.

---

## Input Requirements

### 1. FASTQ Files

- Raw reads from nanopore sequencing (single or multiple `.fastq.gz` files per sample).
- Organized either in subfolders named by barcode or as single files named with the sample ID.

### 2. Barcode Information File (`barcodes.txt`)

A CSV with **no header**, containing:

```{text}
barcode,sample_id,fragments
```

Example:

```{text}
barcode01,SAMPLE001,"Fragment_1 Fragment_2 Fragment_3"
barcode02,SAMPLE002,"FL"
```

### 3. Reference Data

- Full genome reference FASTA with "barcode" region included.
- BED file specifying region coordinates (used to break reference into expected contigs).
- Additional resources in `ReferenceSequences/`:
  - `SIVMac239FullGenome_wBarcode.fas`
  - `SIVregions_wBarcode.bed`
  - `SIVprimers.fasta`
  - `barcode_reference.fas`

---

## Pipeline Steps

1. **Concatenate FASTQ**
   - Merges multiple read files per sample or uses a single file if available.

2. **Quality Filtering**
   - Filters reads using `fastplong` based on quality and expected fragment sizes.

3. **Assembly (Flye)**
   - De novo assembly for each fragment using `flye`.
   - Falls back to using the reference if assembly fails.

4. **Contig Cleanup**
   - Reorders and renames contigs to match reference fragment order.

5. **Polishing**
   - Polishes assemblies using `medaka` (GPU support available).

6. **Read Alignment**
   - Aligns reads back to polished consensus using `minimap2` and `samtools`.

7. **Consensus Trimming**
   - Trims consensus sequences by coverage depth.
   - Special handling for `Fragment_3` in SIV (splits at barcode).

8. **Haplotype Analysis**
   - Counts reads per haplotype using a custom script.

9. **Barcode Counting** *(Optional)*
   - Aligns to barcode reference and counts barcodes using a Perl script.

---

## Configuration Parameters

Defined at the top of the Nextflow script or passed via `-params-file`.

| Parameter | Description | Default |
|----------|-------------|---------|
| `barcodes` | Path to barcode/sample info file | `"barcodes.txt"` |
| `fastq_dir` | Path to directory of FASTQ files | `""` |
| `regions_bed` | BED file of region coordinates | `ReferenceSequences/SIVregions_wBarcode.bed` |
| `reference` | Full genome reference FASTA | `ReferenceSequences/SIVMac239FullGenome_wBarcode.fas` |
| `outdir` | Output directory | `nf-results/` |
| `virus` | Virus name, affects trimming logic | `"SIV"` |
| `min_read_length` | Minimum read length to retain | `-1` (auto) |
| `max_read_length` | Maximum read length to retain | `-1` (auto) |
| `count_barcodes` | Enable barcode count logic | `false` |
| `gpu` | Use GPU for Medaka polishing | `false` |

---

## Usage

### Basic Run

#### Run directly from github

```bash
nextflow rlr-lab/haplotypes_workflow \
  --fastq_dir /path/to/fastqs \
  --barcodes barcodes.txt \
  -work-dir /path/for/work
```

#### If the pipeline is downloaded locally

```bash
nextflow run main.nf \
  --fastq_dir /path/to/fastqs \
  --barcodes barcodes.txt \
  -work-dir /path/for/work
```

### Using a Parameters File

Create a `params.config` file and run:

```bash
nextflow run main.nf -params-file params.config
```

---

## Dependencies

All dependencies are managed through `conda`. The environment is defined in `medaka.yaml`.

This pipeline is designed for execution on a SLURM-based HPC cluster, but can be adapted to other environments.

---

## Output

Results are organized as:

```{text}
nf-results/
  └── SAMPLE_ID/
      └── Fragment_X/
          ├── fastplong/
          │   ├── *_filtered.fastq.gz
          │   ├── fastp_report.*
          ├── flye_out/
          │   └── assembly.fasta
          ├── *_consensus.fasta
          ├── *_aligned.bam(.bai)
          ├── *_consensus_trimmed.fasta
          ├── read_counts.txt
```

There will also be a HTML report generated in the directory the pipeline is launched from. This report contains metrics about workflow execution and contains three main sections: `Summary`, `Resources`, and `Tasks`.

---

## Notes

- The pipeline supports both full-length and tiled fragment designs.
- Barcode alignment is optional and only runs for barcoded fragments.
- For `SIV Fragment_3`, trimmed sequences are further split at a defined barcode delimiter.
- Custom Python and Perl scripts are located in the `scripts/` directory.

---

## Authors & Acknowledgments

Developed by the RLR Lab at Northwestern University.  
Inspired by previous pipelines for tiled amplicon assembly.  
Uses `Nextflow`, `Flye`, `Medaka`, `minimap2`, `samtools`, `fastplong`.

---

## License

This pipeline is released under the GPL-3.0 License.
