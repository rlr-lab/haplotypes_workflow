#!/usr/bin/env python3
import pysam
import numpy as np
import argparse
from collections import defaultdict

def compute_lengths_by_contig(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig_lengths = defaultdict(list)

    for read in bam:
        if read.is_unmapped:
            continue
        contig = bam.get_reference_name(read.reference_id)
        aligned_len = sum(length for op, length in read.cigartuples if op in (0, 2, 3, 7, 8))  # M, D, N, =, X
        contig_lengths[contig].append(aligned_len)

    bam.close()
    return contig_lengths

def filter_bam_by_contig_length(bam_path, output_path, thresholds):
    bam = pysam.AlignmentFile(bam_path, "rb")
    out_bam = pysam.AlignmentFile(output_path, "wb", template=bam)

    for read in bam:
        if read.is_unmapped:
            continue
        contig = bam.get_reference_name(read.reference_id)
        aligned_len = sum(length for op, length in read.cigartuples if op in (0, 2, 3, 7, 8))
        if aligned_len >= thresholds.get(contig, 0):
            out_bam.write(read)

    bam.close()
    out_bam.close()

def main():
    parser = argparse.ArgumentParser(description="Filter BAM reads per contig by dropping shortest X percent.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file after filtering")
    parser.add_argument("--drop-shortest-percent", type=float, default=10,
                        help="Percent of shortest reads to drop (per contig, default: 10)")

    args = parser.parse_args()

    print("Computing read lengths by contig...")
    contig_lengths = compute_lengths_by_contig(args.input_bam)

    thresholds = {
        contig: int(np.percentile(lengths, args.drop_shortest_percent))
        for contig, lengths in contig_lengths.items()
    }

    for contig, cutoff in thresholds.items():
        print(f"Contig {contig}: threshold = {cutoff} bp")

    print("Filtering reads per contig...")
    filter_bam_by_contig_length(args.input_bam, args.output_bam, thresholds)

    print("Done.")

if __name__ == "__main__":
    main()
