#!/usr/bin/env python3
import argparse
import pysam
from collections import defaultdict, Counter

def reads_covering_entire_contig(
    bam_path,
    output_bam=None,
    output_fastq=None,
    min_quality=0,
    min_mapping_quality=0
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    out_bam = None
    out_fastq = open(output_fastq, "w") if output_fastq else None
    contig_read_counts = defaultdict(Counter)

    if output_bam:
        out_bam = pysam.AlignmentFile(output_bam, "wb", template=bam)

    for contig in bam.references:
        print(f'Processing contig {contig}...')
        contig_len = bam.get_reference_length(contig)

        for read in bam.fetch(contig):
            if read.is_unmapped or read.query_sequence is None:
                continue
            if read.mapping_quality < min_mapping_quality:
                continue

            ref_positions = read.get_reference_positions(full_length=True)
            query_sequence = read.query_sequence
            query_qualities = read.query_qualities

            # Filter positions that are in contig bounds and meet quality threshold
            in_bounds = [
                (i, ref_pos)
                for i, ref_pos in enumerate(ref_positions)
                if (
                    ref_pos is not None and
                    0 <= ref_pos < contig_len and
                    query_qualities[i] >= min_quality
                )
            ]

            if not in_bounds:
                continue

            covered_positions = {ref_pos for _, ref_pos in in_bounds}

            # Check that all positions in the contig are covered
            if not all(pos in covered_positions for pos in range(contig_len)):
                continue  # Read no longer fully covers the contig

            trimmed_seq = ''.join(query_sequence[i] for i, _ in in_bounds)
            trimmed_qual_str = ''.join(chr(query_qualities[i] + 33) for i, _ in in_bounds)

            contig_read_counts[contig][trimmed_seq] += 1

            if out_fastq:
                out_fastq.write(f"@{read.query_name}_{contig}\n")
                out_fastq.write(trimmed_seq + "\n+\n" + trimmed_qual_str + "\n")

            if out_bam:
                out_bam.write(read)

    bam.close()
    if out_bam:
        out_bam.close()
    if out_fastq:
        out_fastq.close()

    return contig_read_counts


def write_counts_to_file(counts, output_path):
    with open(output_path, 'w') as f:
        for contig, seq_counts in counts.items():
            for seq, count in seq_counts.items():
                f.write(f"{contig}\t{seq}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Filter and count reads that cover entire contigs, with trimming and quality thresholds."
    )
    parser.add_argument("input_bam", help="Input sorted BAM file")
    parser.add_argument("--output-bam", help="Write filtered reads to a new BAM file")
    parser.add_argument("--output-fastq", help="Write trimmed sequences to a FASTQ file")
    parser.add_argument("--count-output", default="read_counts.txt", help="Write read count table to this file (default = read_counts.txt)")
    parser.add_argument("--min-quality", type=int, default=0, help="Minimum base quality (Phred score, default = 0)")
    parser.add_argument("--min-mapping-quality", type=int, default=20, help="Minimum mapping quality (MAPQ, default = 20)")

    args = parser.parse_args()

    counts = reads_covering_entire_contig(
        bam_path=args.input_bam,
        output_bam=args.output_bam,
        output_fastq=args.output_fastq,
        min_quality=args.min_quality,
        min_mapping_quality=args.min_mapping_quality
    )

    if args.count_output:
        write_counts_to_file(counts, args.count_output)


if __name__ == "__main__":
    main()

