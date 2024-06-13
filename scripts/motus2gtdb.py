#!/usr/bin/env python

import argparse
import os.path
import sys

##### This script transforms mOTUs results from NCBI taxonomy to GTDB
##### Alejandra Escobar, EMBL-EBI
##### May 3, 2024

def motus_parser( motus_table ):
    taxonomy_dir = {}
    with open(motus_table, 'r') as input_file:
        next( input_file )
        for line in input_file:
            motu, consensus_taxonomy, count = line.rstrip().split('\t')
            taxonomy_dir[motu] = count
    return( taxonomy_dir )


def gtdb_mapping( taxonomy_dir, mapping_table, sample_name ):
    gtdb_lineages = {}
    with open(mapping_table, 'r') as input_file:
        for line in input_file:
            line_l = line.rstrip().split('\t')
            motu_id = line_l.pop(0)
            if motu_id in taxonomy_dir:
                clean_lineage = [item for item in line_l if 'Not_annotated' not in item and 'Incongruent' not in item]
                lineage = ';'.join(clean_lineage)
                abundance  =  int(taxonomy_dir[motu_id])
                if lineage in gtdb_lineages:
                    new_abundance = gtdb_lineages[lineage] + abundance
                    gtdb_lineages[lineage] = new_abundance
                else:
                    gtdb_lineages[lineage] = abundance


    with open(sample_name + '_gtdb.tsv', 'w') as output_file:
        output_file.write('\t'.join(['lineage', sample_name]) + '\n' )
        for tax_name in gtdb_lineages:
            output_file.write('\t'.join([tax_name, str(gtdb_lineages[tax_name])]) + '\n')



def main():
    parser = argparse.ArgumentParser(
            description="This script transforms mOTUs results from NCBI taxonomy to GTDB"
    )
    parser.add_argument(
        "--input",
        type=str,
        help="mOTUs result table with NCBI taxonomy (*_merged.fastq.tsv)",
        required=True,
    )
    parser.add_argument(
        "--mapping",
        type=str,
        help="Mapping file with GTDB taxonomy (mOTUs_3.0.0_GTDB_tax.tsv)",
        required=True,
    )
    parser.add_argument(
       "--sample",
        type=str,
        help="Sample name prefix to use in the header line",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    taxonomy_dir = motus_parser( args.input )
    gtdb_mapping( taxonomy_dir, args.mapping, args.sample )

if __name__ == "__main__":
    main()

