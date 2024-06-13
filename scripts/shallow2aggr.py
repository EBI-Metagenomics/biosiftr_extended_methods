#!/usr/bin/env python

import argparse
import os.path
import sys

##### This script clean and aggregate lineages from shallo_mapping results
##### Alejandra Escobar, EMBL-EBI
##### May 8, 2024


def tax_parser( taxonomy_table ):
    aggregated = {}
    with open(taxonomy_table, 'r') as input_file:
        header = input_file.readline().strip().replace('#OTU ID','feature')
        for line in input_file:
            line_l = line.rstrip().split('\t')
            lineage = line_l.pop(0)
            lineage = lineage.split(';')
            lineage.pop(-1)
            filtered_lineage = [item for item in lineage if len(item.strip()) != 3]
            filtered_lineage = ';'.join(filtered_lineage)

            current_values = [float(item) for item in line_l]

            if filtered_lineage in aggregated:
                saved_values = aggregated[filtered_lineage]
                added_values = [x + y for x, y in zip(saved_values, current_values)]
                aggregated[filtered_lineage] = added_values
            else:
                aggregated[filtered_lineage] = current_values
            
    # Writing output
    with open('aggregated.tsv', 'w') as output_file:
        output_file.write(header + '\n')
        for lineage in aggregated:
            values = [str(item) for item in aggregated[lineage]]
            output_file.write( lineage + '\t' + '\t'.join(values)+ '\n')



def main():
    parser = argparse.ArgumentParser(
            description="This script clean and aggregate lineages from shallow_mapping results"
    )
    parser.add_argument(
        "--taxonomy_table",
        type=str,
        help="Shallow-mapping integrated taxonomy table ([sm|bwa]_taxo_matrix.tsv)",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    tax_parser( args.taxonomy_table )


if __name__ == "__main__":
    main()

