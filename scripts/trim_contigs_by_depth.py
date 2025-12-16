#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_coverage_regions(bam, contig_id, min_depth):
    """Returns the start and end of the high-coverage region for a contig."""
    try:
        cmd = ["samtools", "depth", "-a", "-r", contig_id, bam]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except:
        return None, None
    
    positions = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        chrom, pos, depth = line.split("\t")
        if int(depth) >= min_depth:
            positions.append(int(pos))

    if not positions:
        return None, None

    return min(positions), max(positions)

def trim_contigs_by_depth(fasta, bam, output_fasta, min_depth):
    records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    trimmed_records = []

    for contig_id, record in records.items():
        print(f"Processing {contig_id}...")
        start, end = get_coverage_regions(bam, contig_id, min_depth)
        if start is None or end is None:
            print(f"Skipping {contig_id} — no region meets depth threshold {min_depth}")
            continue
        trimmed_seq = record.seq[start-1:end]  # convert to 0-based
        new_record = SeqRecord(
            trimmed_seq,
            id=record.id,
            description=f"trimmed {start}-{end} (min depth {min_depth})"
        )
        trimmed_records.append(new_record)

    SeqIO.write(trimmed_records, output_fasta, "fasta")
    print(f"Trimmed contigs written to {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Trim contigs to regions with high coverage based on BAM depth.")
    parser.add_argument("-f", "--fasta", required=True, help="Input contig FASTA file")
    parser.add_argument("-b", "--bam", required=True, help="Sorted and indexed BAM file")
    parser.add_argument("-o", "--output", default="trimmed_contigs.fasta", help="Output trimmed FASTA")
    parser.add_argument("-d", "--depth", type=float, default=1000, help="Minimum depth threshold (default: 1000)")
    args = parser.parse_args()

    trim_contigs_by_depth(args.fasta, args.bam, args.output, args.depth)

if __name__ == "__main__":
    main()
