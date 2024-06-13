#!/usr/bin/env python

import argparse
import os.path
import sys

##### This script transforms a counts table into an absence-presence table
##### Alejandra Escobar, EMBL-EBI
##### May 16, 2024


def matrix_parser( matrix, output ):
    with open(matrix, 'r') as input_file, open(output, 'w') as output_file:
        header = input_file.readline().strip()
        output_file.write(header + '\n')
        for line in input_file:
            line_l = line.rstrip().split('\t')
            feature = line_l.pop(0)
            values = [float(item) for item in line_l]
            values_t = [1 if x > 0 else 0 for x in values]
            values_t = [str(item) for item in values_t]
            output_file.write(feature + '\t' + '\t'.join(values_t) + '\n')



def main():
    parser = argparse.ArgumentParser(
            description="This script transforms a counts table into an absence-presence table"
    )
    parser.add_argument(
        "--matrix",
        type=str,
        help="Input file to be transformed into absence-presence",
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
    matrix_parser( args.matrix, args.output )


if __name__ == "__main__":
    main()

