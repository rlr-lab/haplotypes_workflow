#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
import os
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ContigHit = namedtuple('ContigHit', ['contig', 'strand', 'ref_start'])

def run_minimap2(reference, assembly, paf_out):
    print("Running minimap2...")
    subprocess.run([
        "minimap2", "-x", "asm5", reference, assembly
    ], stdout=open(paf_out, 'w'), check=True)

def parse_paf(paf_file):
    hits = {}
    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            contig = parts[0]
            strand = parts[4]
            ref_start = int(parts[7])
            if contig not in hits or ref_start < hits[contig].ref_start:
                hits[contig] = ContigHit(contig, strand, ref_start)
    return sorted(hits.values(), key=lambda x: x.ref_start)

def reorder_and_reorient_contigs(assembly_fasta, hits, output_fasta, concat=False, linker=""):
    contig_dict = SeqIO.to_dict(SeqIO.parse(assembly_fasta, "fasta"))

    if concat:
        assembled_seq = ""
        contig_ids = []
        for hit in hits:
            record = contig_dict[hit.contig]
            seq = record.seq.reverse_complement() if hit.strand == "-" else record.seq
            assembled_seq += str(seq) + linker
            contig_ids.append(hit.contig)

        #if assembled_seq.endswith(linker):
        #    assembled_seq = assembled_seq[:-len(linker)]

        single_record = SeqRecord(
            seq=Seq(assembled_seq),
            id="assembled_sequence",
            description=f"Joined from: {' '.join(contig_ids)}"
        )
        SeqIO.write([single_record], output_fasta, "fasta")
    else:
        ordered_records = []
        for hit in hits:
            record = contig_dict[hit.contig]
            if hit.strand == "-":
                record.seq = record.seq.reverse_complement()
                record.description += " (reverse complemented)"
            ordered_records.append(record)
        SeqIO.write(ordered_records, output_fasta, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Reorder and reorient contigs against a reference using minimap2.")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome (.fasta)")
    parser.add_argument("-a", "--assembly", required=True, help="De novo contigs (.fasta)")
    parser.add_argument("-o", "--output", default="contigs.sorted.fasta", help="Output FASTA file")
    parser.add_argument("--concat", action="store_true", help="Concatenate contigs into one sequence with Ns between")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        paf_path = os.path.join(tmpdir, "alignment.paf")
        run_minimap2(args.reference, args.assembly, paf_path)
        hits = parse_paf(paf_path)
        reorder_and_reorient_contigs(args.assembly, hits, args.output, concat=args.concat)

    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()

