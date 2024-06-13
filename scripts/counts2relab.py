#!/usr/bin/env python

import argparse
import os.path
import sys

##### This script transforms a counts table into relative abundance and removes singletones and doubletons
##### Alejandra Escobar, EMBL-EBI
##### May 8, 2024


def count_parser( count_table, output ):
    lineage_values = {}
    singleton_counter = 0
    with open(count_table, 'r') as input_file:
        header = input_file.readline().strip()
        len_total_sums = len(header) - 1
        total_sums = [0] * len_total_sums
        for line in input_file:
            line_l = line.rstrip().split('\t')
            lineage = line_l.pop(0).replace(' ','_')
            values = [float(item) for item in line_l]
            val_sum = sum(values)
            if val_sum > 2:
                total_sums = [x + y for x, y in zip(total_sums, values)]
                lineage_values[lineage] = values
            else:
                singleton_counter += 1
    print('Discarded singletons and doubletons: '+str(singleton_counter))

    # Calculating relative abundace and writing output
    with open(output, 'w') as output_file:
        output_file.write(header + '\n')
        for lineage in lineage_values:
            rel_ab = [x / y for x, y in zip(lineage_values[lineage], total_sums)]
            vals_str = [str(item) for item in rel_ab]
            output_file.write( lineage + '\t' + '\t'.join(vals_str)+ '\n')


def main():
    parser = argparse.ArgumentParser(
            description="This script transforms a counts matrix into relative abundance and removes singletones and doubletons"
    )
    parser.add_argument(
        "--count_table",
        type=str,
        help="Count table to be transformed into relative abundance",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output file name",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    count_parser( args.count_table, args.output )


if __name__ == "__main__":
    main()

