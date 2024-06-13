#!/usr/bin/env python

import argparse
import os.path
import numpy as np
import sys

##### This script computes precision, recall, and F1 score from a counts matrix
##### Alejandra Escobar, EMBL-EBI
##### Jan 17, 2024


def matrix_parser( matrix_file):
    bwa_dict = {}
    sm_dic = {}
    mock_dic = {}
    all_features = []
    with open(matrix_file, 'r') as input_file:
        samples_list = input_file.readline().strip().split('\t')
        samples_list.pop(0)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            feature = l_line.pop(0)
            all_features.append(feature)
            index = 0
            for element in l_line:
                sample_name = samples_list[index]
                sample_group = sample_name.split('_')[0]
                index+=1
                if sample_group == 'GT':
                    short_sample_name = sample_name.replace('GT_','')
                    composite_key = (short_sample_name, feature)
                    mock_dic[composite_key] = float(element)
                elif sample_group == 'sourmash':
                    short_sample_name = sample_name.replace('sourmash_','')
                    composite_key = (short_sample_name, feature)
                    sm_dic[composite_key] = float(element)
                elif sample_group == 'bwamem2':
                    short_sample_name = sample_name.replace('bwamem2_','')
                    composite_key = (short_sample_name, feature)
                    bwa_dict[composite_key] = float(element)

    return( bwa_dict, sm_dic, mock_dic, all_features )


def confusion_counter( bwa_dict, sm_dic, mock_dic, all_features ):
    performance_dict = {}
    for num in range(1,21):
        dataset = str(num)
        performance_dict[(dataset, 'bwamem2', 'TP')] = 0
        performance_dict[(dataset, 'bwamem2', 'FP')] = 0
        performance_dict[(dataset, 'bwamem2', 'FN')] = 0
        performance_dict[(dataset, 'sourmash', 'TP')] = 0
        performance_dict[(dataset, 'sourmash', 'FP')] = 0
        performance_dict[(dataset, 'sourmash', 'FN')] = 0
        for feature in all_features:
            composite_key = (dataset, feature)
            ref_value = mock_dic[composite_key]

            # Parsing bwa predictions
            obs_value = bwa_dict[composite_key]
            if all([ ref_value>0 , obs_value>0 ]):
                performance_dict[(dataset, 'bwamem2', 'TP')]+=1
            if all([ ref_value==0 , obs_value>0 ]):
                performance_dict[(dataset, 'bwamem2', 'FP')]+=1
            if all([ ref_value>0 , obs_value==0 ]):
                performance_dict[(dataset, 'bwamem2', 'FN')]+=1

            # Parsing sourmash predictions
            obs_value = sm_dic[composite_key]
            if all([ ref_value>0 , obs_value>0 ]):
                performance_dict[(dataset, 'sourmash', 'TP')]+=1
            if all([ ref_value==0 , obs_value>0 ]):
                performance_dict[(dataset, 'sourmash', 'FP')]+=1
            if all([ ref_value>0 , obs_value==0 ]):
                performance_dict[(dataset, 'sourmash', 'FN')]+=1

    return(performance_dict)


def table_writer( performance_dict, output_name ):
    tools = ['bwamem2', 'sourmash']
    with open( output_name, 'w') as output_file:
        output_file.write("\t".join(['dataset','predictor','metric','value'])+'\n')
        for num in range(1,21):
            for tool in tools:
                (precision, recall, f1_score) = metrics_calc(performance_dict, str(num), tool)
                dataset = 'synth_'+str(num)
                output_file.write("\t".join([dataset,tool,'precision',str(precision)])+'\n')
                output_file.write("\t".join([dataset,tool,'recall',str(recall)])+'\n')
                output_file.write("\t".join([dataset,tool,'f1_score',str(f1_score)])+'\n')



def metrics_calc(performance_dict, dataset, tool):
    TP = performance_dict[(dataset, tool, 'TP')]
    FP = performance_dict[(dataset, tool, 'FP')]
    FN = performance_dict[(dataset, tool, 'FN')]
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    f1_score = 2 * ( precision * recall ) / ( precision + recall )
    return( precision, recall, f1_score )



def main():
    parser = argparse.ArgumentParser(
            description="This script computes precision, recall, and F1 score from a counts matrix"
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Counts matrix",
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
    (bwa_dict, sm_dic, mock_dic, all_features) = matrix_parser( args.input )
    (performance_dict) = confusion_counter( bwa_dict, sm_dic, mock_dic, all_features )
    table_writer( performance_dict, args.output )

if __name__ == "__main__":
    main()

