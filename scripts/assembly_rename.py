#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from Bio import SeqIO
import sys
import argparse
import os

##### This script filter a fasta file by a given sequence size and rename the contigs
##### Alejandra Escobar, EMBL-EBI
##### January 18th, 2019

parser = argparse.ArgumentParser(
        description='This script filter a fasta file by a given sequence size and rename the contigs')
parser.add_argument('assembly.fasta', type=str, help='Input fasta file')
parser.add_argument('prefix', type=str, help='Prefix to be added to the contig name')
args = parser.parse_args()

### Reading the fasta file
input_file=sys.argv[1]
prefix=sys.argv[2]
output_file='renamed_'+prefix+'.fasta'

counter=0
with open(output_file, 'w') as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        my_chain=str(record.seq).upper()
        counter+=1
        output_handle.write(">"+prefix+"_"+str(counter)+'\n')
        output_handle.write(my_chain+'\n')


