#!/usr/bin/env python

import argparse
import os.path
import sys
import re

##### This script cross ASVs ids with their corresponding taxonomic labels
##### Alejandra Escobar, EMBL-EBI
##### May 3, 2024


def asv_parser( asv_table ):
    asv_dir = {}
    with open(asv_table, 'r') as input_file:
        next( input_file )
        header = input_file.readline().strip().replace('#OTU ID','feature')
        for line in input_file:
            line_l = line.rstrip().split('\t')
            otu_id = line_l.pop(0)
            values = [float(item) for item in line_l]
            asv_dir[otu_id] = values
    return( asv_dir, header )


def tax_parser( names_table ):
    pattern_1 = r'\w+_\d+$'
    pattern_2 = r'_([A-Z])_([0-9]+) ([a-z]+)'
    tax_dir = {}
    with open(names_table, 'r') as input_file:
        next( input_file )
        for line in input_file:
            asv, lineage, conf = line.rstrip().split('\t')
            lineage = lineage.replace('; ',';').split(';')
            filtered_lineage = [item for item in lineage if len(item.strip()) != 3]

            new_lineage = []
            for rank in filtered_lineage:
                ## Removing monophyletic identifiers in intermediate ranks
                if re.search(pattern_1, rank):
                    rank = rank.split('_')
                    rank.pop(-1)
                    rank = '_'.join(rank)
                    new_lineage.append(rank)
                else:
                    ## Removing monophyletic identifiers in genus
                    if all([ rank.startswith('s__'), re.search(pattern_2, rank) ]):
                        tagged_name, species = rank.split(' ')
                        tagged_name = tagged_name.split('_')
                        tagged_name.pop(-1)
                        clean_rank = '_'.join(tagged_name)+' '+species
                        new_lineage.append(clean_rank)
                    else:
                        new_lineage.append(rank)

            new_lineage = ';'.join(new_lineage)
            filtered_lineage = ';'.join(new_lineage).replace(' ','_')
            tax_dir[asv] = new_lineage

    return( tax_dir )


def output_writer( asv_dir, tax_dir, header ):
    # Aggregating counts per taxonomic label
    aggregated = {}
    for asv_id in asv_dir:
        if asv_id in tax_dir:
            current_values = asv_dir[asv_id]
            lineage = tax_dir[asv_id]
            if lineage in aggregated:
                saved_values = aggregated[lineage]
                added_values = [x + y for x, y in zip(saved_values, current_values)]
                aggregated[lineage] = added_values
            else:
                aggregated[lineage] = current_values
        else:
            print('No lineage for '+asv_id+'\n')
            
    # Writing output
    with open('lineages.tsv', 'w') as output_file:
        output_file.write(header + '\n')
        for lineage in aggregated:
            values = [str(item) for item in aggregated[lineage]]
            output_file.write( lineage + '\t' + '\t'.join(values)+ '\n')



def main():
    parser = argparse.ArgumentParser(
            description="This script cross ASVs ids with their corresponding taxonomic labels"
    )
    parser.add_argument(
        "--asv_table",
        type=str,
        help="ASV abundance table with GG2 taxonomy (ASV-table_gg2.tsv)",
        required=True,
    )
    parser.add_argument(
        "--names_table",
        type=str,
        help="Mapping file from ASVs to lineages (taxonomy.tsv)",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( asv_dir, header ) = asv_parser( args.asv_table )
    tax_dir = tax_parser( args.names_table )
    output_writer( asv_dir, tax_dir, header )

if __name__ == "__main__":
    main()

