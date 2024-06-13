#!/usr/bin/env python

import argparse
import os.path
import sys
import re

##### This script aggregates a taxonomic matrix to all 7 taxonomic levels
##### Alejandra Escobar, EMBL-EBI
##### May 30, 2024


def tax_parser( tax_table ):
    tax_dict = {}
    with open(tax_table, 'r') as input_file:
        header = input_file.readline().strip()
        for line in input_file:
            line_l = line.rstrip().split('\t')
            lineage = line_l.pop(0)
            values = [float(item) for item in line_l]
            tax_dict[lineage] = values
    return( tax_dict, header )


def add_ranks( composite_key, aggregated_dict, current_values ):
    if composite_key in aggregated_dict:
        saved_values = aggregated_dict[composite_key]
        added_values = [x + y for x, y in zip(saved_values, current_values)]
        aggregated_dict[composite_key] = added_values
    else:
        aggregated_dict[composite_key] = current_values

    return aggregated_dict


def lineage_breaker( tax_dict ):
    # Aggregating counts per taxonomic label
    aggregated_dict = {}
    tax_levels = {
        0 : 'domain',
        1 : 'phylum',
        2 : 'class',
        3 : 'order',
        4 : 'family',
        5 : 'genus',
        6 : 'species'
    }
    for lineage in tax_dict:
        current_values = tax_dict[lineage]
        ranks = lineage.split(';')

        # Aggregating at each taxonomic level
        for index in tax_levels:
            if len(ranks) >= index+1:
                composite_key = ( tax_levels[index], ranks[index] )
                aggregated_dict = add_ranks( composite_key, aggregated_dict, current_values )

    return(aggregated_dict)


def outputs_writer( aggregated_dict, header, out_prefix ):
    ranks_list = [
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species'
    ]

    # Writing output files
    for tax_level in ranks_list:
        output_name = out_prefix + '_' + tax_level + '.tsv'
        with open(output_name, 'w') as output_file:
            output_file.write(header + '\n')
            for composite_key in aggregated_dict:
                if composite_key[0] == tax_level:
                    rounded_numbers = [round(num, 4) for num in aggregated_dict[composite_key]]
                    values = [str(item) for item in rounded_numbers]
                    output_file.write( composite_key[1] + '\t' + '\t'.join(values)+ '\n')



def main():
    parser = argparse.ArgumentParser(
            description="This script aggregates a taxonomic matrix to all 7 taxonomic levels"
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Taxonomic tsv table. Lineage in the first colum",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output prefix. A rank level suffux will be added to each output. Example: relab_table",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( tax_dict, header ) = tax_parser( args.input )
    aggregated_dict = lineage_breaker( tax_dict )
    outputs_writer( aggregated_dict, header, args.output )

if __name__ == "__main__":
    main()

