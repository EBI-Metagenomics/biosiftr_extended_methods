#!/usr/bin/env python

import sys
import gzip
import argparse
from Bio import SeqIO

##### Alejandra Escobar, EMBL-EBI/Sanger
##### Aug 9th, 2021


### Input files
parser = argparse.ArgumentParser(
        description='This script transform fastq to fasta modifying the header')
parser.add_argument('input.fastq', type=str, help='fastq file')
parser.add_argument('read_prefix', type=str, help='reads prefix')
args = parser.parse_args()

input_file = sys.argv[1]
read_prefix = sys.argv[2] 
counter=0

with gzip.open(input_file, "rt") as handle, open(read_prefix+'.fasta','w') as output:
    for record in SeqIO.parse(handle, "fastq"):
        counter+=1
        actual_seq=str(record.seq).upper()
        output.write('>'+read_prefix+'_'+str(counter)+'\n')
        output.write(str(record.seq).upper()+'\n')

