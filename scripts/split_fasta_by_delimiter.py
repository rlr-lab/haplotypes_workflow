#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# IUPAC ambiguity codes mapped to regex patterns
IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]"
}

def iupac_to_regex(delimiter):
    """Convert an IUPAC DNA sequence to a regex pattern."""
    return "".join(IUPAC_CODES.get(base.upper(), base) for base in delimiter)

def split_fasta_by_delimiter(input_fasta, delimiter_seq, output_fasta):
    pattern = re.compile(iupac_to_regex(delimiter_seq), re.IGNORECASE)
    output_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        full_seq = str(record.seq)
        match = pattern.search(full_seq)

        if not match:
            print(f"No match found in {record.id}. Skipping.")
            continue

        start, end = match.span()
        seq1 = full_seq[:start]
        seq2 = full_seq[end:]

        rec1 = SeqRecord(
            seq=Seq(seq1),
            id=f"{record.id}_part1",
            description="first part before delimiter"
        )
        rec2 = SeqRecord(
            seq=Seq(seq2),
            id=f"{record.id}_part2",
            description="second part after delimiter"
        )

        output_records.extend([rec1, rec2])
        break  # stop after first sequence is split

    if output_records:
        SeqIO.write(output_records, output_fasta, "fasta")
        print(f"Wrote 2 split contigs to: {output_fasta}")
    else:
        print("No valid sequences were split and written.")

def main():
    parser = argparse.ArgumentParser(description="Split a FASTA sequence using an ambiguous DNA delimiter.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-d", "--delimiter", required=True, help="Ambiguous DNA sequence delimiter (e.g. ACGNNNTTA)")
    parser.add_argument("-o", "--output", default="split_contigs.fasta", help="Output FASTA file containing both split parts")
    args = parser.parse_args()

    split_fasta_by_delimiter(args.input, args.delimiter, args.output)

if __name__ == "__main__":
    main()
