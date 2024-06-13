#!/usr/bin/env python

import argparse
import os.path
import gzip
import sys
from Bio import SeqIO

##### This script extracts KOs annotation from emapper results using abe eval threshold = 1e-10
##### Alejandra Escobar, EMBL-EBI
##### May 16, 2024


def eggnog_parser( eggnog_file, sample_name ):
    kos_list = {}
    with open(eggnog_file, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            evalue = float(l_line[2])
            if evalue <= 1e-10:
                kegg_ko = l_line[11].replace('ko:','')
                if kegg_ko != '-':
                    for ko in kegg_ko.split(','):
                        if ko in kos_list:
                            kos_list[ko]+=1
                        else:
                            kos_list[ko]=1

    with open(sample_name+'_kos.txt', 'w') as output_file:
        output_file.write('ko\t'+sample_name+'\n')
        for ko in kos_list:
            output_file.write(ko+'\t'+str(kos_list[ko])+'\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script extracts KOs annotation from emapper results using abe eval threshold = 1e-10"
    )
    parser.add_argument(
        "--eggnog",
        type=str,
        help="EggNOG results path",
        required=True,
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Sample name to be used in the header and to name the output file. Default = sample",
        required=False,
    )
    args = parser.parse_args()

    if args.sample:
        output_name = args.sample
    else:
        output_name = 'sample'

    eggnog_parser( args.eggnog, output_name )


if __name__ == "__main__":
    main()

