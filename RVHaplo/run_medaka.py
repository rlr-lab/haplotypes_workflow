import os
import sys
import numpy as np
from collections import Counter
import pickle

file_path = sys.argv[1]
file_prefix = sys.argv[2]
flag = sys.argv[3]

if flag == '1':
    os.system(f"rm -rf {file_path}/medaka/medaka0")
    os.system(f"mkdir {file_path}/medaka/medaka0")
    # os.system(f'medaka_consensus -i {file_path}/medaka/fastx/cluster.fastq -d {file_path}/medaka/fastx/consensus.fasta -o {file_path}/medaka/medaka0')
    os.system(f'mini_align -i {file_path}/medaka/fastx/cluster.fastq -r {file_path}/medaka/fastx/consensus.fasta -m -t 4 -p {file_path}/medaka/reads')
    os.system(f'medaka inference {file_path}/medaka/reads.bam {file_path}/medaka/contig.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0')
    os.system(f'medaka sequence {file_path}/medaka/contig.hdf {file_path}/medaka/fastx/consensus.fasta {file_path}/medaka/medaka0/consensus.fasta --threads 4')
    os.system(f'cp {file_path}/medaka/medaka0/consensus.fasta {file_path}/{file_prefix}_haplotypes.fasta')
else:
    f = open(f'{file_path}/{file_prefix}_clusters.pickle','rb')
    temp = Counter(pickle.load(f)['haplo_fre'])
    num = len(temp)
    f.close()
    for i in range(num):
        os.system(f"rm -rf {file_path}/medaka/medaka_{i}")
        os.system(f'mkdir {file_path}/medaka/medaka_{i}')
        # os.system(f'medaka_consensus -i {file_path}/medaka/fastx/cluster_{i}.fastq -d {file_path}/medaka/fastx/consensus_{i}.fasta -o {file_path}/medaka/medaka_{i}')
        os.system(f'mini_align -i {file_path}/medaka/fastx/cluster_{i}.fastq -r {file_path}/medaka/fastx/consensus_{i}.fasta -m -t 4 -p {file_path}/medaka/medaka_{i}/reads')
        os.system(f'medaka inference {file_path}/medaka/medaka_{i}/reads.bam {file_path}/medaka/medaka_{i}/contig.hdf --threads 2 --model r1041_e82_400bps_hac_v4.3.0')
        os.system(f'medaka sequence {file_path}/medaka/medaka_{i}/contig.hdf {file_path}/medaka/fastx/consensus_{i}.fasta {file_path}/medaka/medaka_{i}/consensus.fasta --threads 4')
        os.system(f"cat {file_path}/medaka/medaka_{0}/consensus.fasta > {file_path}/{file_prefix}_haplotypes.fasta")
    for i in range(1,num):
        os.system(f"cat {file_path}/medaka/medaka_{i}/consensus.fasta >> {file_path}/{file_prefix}_haplotypes.fasta")

exit()
