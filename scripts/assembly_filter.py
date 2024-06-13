#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse

##### This script filter a fasta file by a given sequence size
##### Alejandra Escobar, EMBL-EBI/Sanger
##### January 18th, 2019

parser = argparse.ArgumentParser(
        description='This script filter a fasta file by a given sequence size')
parser.add_argument('assembly.fasta', type=str, help='Input fasta file')
parser.add_argument('seq_length', type=int, help='Contigs of sequence length < seq_length will be discarted')
args = parser.parse_args()

### Reading the fasta file
input_file=sys.argv[1]
threshold=int(sys.argv[2])

for record in SeqIO.parse(input_file, "fasta"):
    my_chain=str(record.seq).upper()
    if len(my_chain)>=threshold:
        print(">"+str(record.id))
        print(my_chain)


