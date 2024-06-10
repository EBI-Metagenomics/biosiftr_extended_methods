#!/usr/bin/env python

import argparse
import os.path
import pysam
import sys
from Bio import SeqIO

##### This script compute precision and recall of predictions based on taxonomic nodes comparison
##### Alejandra Escobar, EMBL-EBI
##### Aug 3, 2023


def metadata_parser( catalogue_metadata ):
    catalogue_lineages = {}
    with open(catalogue_metadata, 'r') as input_file:
        next(input_file)
        for line in input_file:
            l_line = line.rstrip().split('\t')
            genome = l_line[0]
            lineage = l_line[14].replace(' ','_')
            catalogue_lineages[genome]=lineage
    return( catalogue_lineages )


def mock_parser( mock_data, expected_covs ):
    mock_nodes = []
    with open(expected_covs, 'r') as input_file:
        next(input_file)
        for line in input_file:
            (
                genome,
                dataset,
                total_reads,
                genome_size,
                assem_abundance,
                genome_coverage,
                mock_lineage) = line.rstrip().split('\t')
            if dataset == mock_data:
                current_nodes = mock_lineage.split(';')
                mock_nodes = mock_nodes + current_nodes
    mock_nodes = set(mock_nodes)
    return( mock_nodes )

def bwa_parser( bwa_table, mapping_type, catalogue_lineages, bwa_nodes ):
    bwa_nodes_pre = {}
    cov_threshold = [
        0,
        0.0001,
        0.001,
        0.01,
        0.1,
        1.0,
        10,
    ]
    with open(bwa_table, 'r') as input_file:
        next(input_file)
        for line in input_file:
            (
                reference_id,
                reference_len,
                reads_count,
                base_cov,
                rel_abun) = line.rstrip().split('\t')

            for thres in cov_threshold:
                if float(base_cov) >= thres:
                    lineage = catalogue_lineages[reference_id].split(';')
                    index = -1
                    index_list = []
                    for rank in lineage:
                        index += 1
                        if len(rank) == 3:
                            index_list.append(index)
                    if len(index_list) > 0:
                        sorted_index = sorted(index_list, reverse=True)
                        for index in sorted_index:
                            lineage.pop(index)

                    composite_key = ( mapping_type, thres )
                    if composite_key in bwa_nodes_pre:
                        bwa_nodes_pre[composite_key] = bwa_nodes_pre[composite_key] + lineage
                    else:
                        bwa_nodes_pre[composite_key] = lineage

    for tuple_key in bwa_nodes_pre:
        clean_nodes = set(bwa_nodes_pre[tuple_key])
        bwa_nodes[tuple_key] = clean_nodes

    return( bwa_nodes )


def sm_parser( sm_table, kmer, sm_nodes ):
    sm_nodes_pre = {}
    with open(sm_table, 'r') as input_file:
        next(input_file)
        for line in input_file:
            comp_lineage, rel_abun = line.rstrip().split('\t')
            lineage = comp_lineage.split(';')
            genome_id = lineage.pop(-1)
            index = -1
            index_list = []
            for rank in lineage:
                index += 1
                if len(rank) == 3:
                    index_list.append(index)
            if len(index_list) > 0:
                sorted_index = sorted(index_list, reverse=True)
                for index in sorted_index:
                    lineage.pop(index)

            if kmer in sm_nodes_pre:
                sm_nodes_pre[kmer] = sm_nodes_pre[kmer] + lineage
            else:
                sm_nodes_pre[kmer] = lineage
    for tuple_key in sm_nodes_pre:
        clean_nodes = set(sm_nodes_pre[tuple_key])
        sm_nodes[tuple_key] = clean_nodes
    return( sm_nodes )

def bwa_metrics_writer( dataset, mock_data, mock_nodes, bwa_nodes ):
    cov_threshold = [
        0,
        0.0001,
        0.001,
        0.01,
        0.1,
        1.0,
        10,
    ]
    bwa_types = [
        'all',
        'best',
        'unique'
    ]
    
    with open('bwa_'+dataset+'_'+mock_data+'_metrics.txt', 'w') as file_out:
        file_out.write('\t'.join([
            'dataset',
            'synth_comm',
            'mapping_type',
            'threshold',
            'metric',
            'value'
        ]) + '\n' )

        gt = mock_nodes
        for bwa_type in bwa_types:
            for threshold in cov_threshold:
                composite_key = ( bwa_type, threshold )
                if composite_key in bwa_nodes:
                    bwa_pred = bwa_nodes[composite_key]
                    ( bwa_recall, bwa_precision, bwa_f1 ) = metrics_calc( gt, bwa_pred )    
                    file_out.write("\t".join([
                        dataset,
                        mock_data,
                        bwa_type,
                        str(threshold),
                        'recall',
                        str(bwa_recall)
                    ])+'\n')
                    file_out.write("\t".join([
                        dataset,
                        mock_data,
                        bwa_type,
                        str(threshold),
                        'precision',
                        str(bwa_precision)
                    ])+'\n')
                    file_out.write("\t".join([
                        dataset,
                        mock_data,
                        bwa_type,
                        str(threshold),
                        'f1_score',
                        str(bwa_f1)
                    ])+'\n')


def sm_metrics_writer( dataset, mock_data, mock_nodes, sm_nodes ):
    sm_kmers = [
        'kmer_11',
        'kmer_21',
        'kmer_31',
        'kmer_51'
    ]

    with open('sm_'+dataset+'_'+mock_data+'_metrics.txt', 'w') as file_out:
        file_out.write('\t'.join([
            'dataset',
            'mock',
            'kmer',
            'metric',
            'value'
        ]) + '\n' )

        gt = mock_nodes
        for sm_kmer in sm_kmers:
            sm_pred = sm_nodes[sm_kmer]
            ( sm_recall, sm_precision, sm_f1 ) = metrics_calc( gt, sm_pred )
            file_out.write("\t".join([
                dataset,
                mock_data,
                sm_kmer,
                'recall',
                str(sm_recall)])+'\n')
            file_out.write("\t".join([
                dataset,
                mock_data,
                sm_kmer,
                'precision',
                str(sm_precision)])+'\n')
            file_out.write("\t".join([
                dataset,
                mock_data,
                sm_kmer,
                'f1_score',
                str(sm_f1)])+'\n')


def metrics_calc( GT, PR ):
    recall = float(len( GT & PR )) / float(len( GT ))
    precision = float(len( GT & PR )) / float(len( PR ))
    f1_score = float(2)*(precision * recall) / (precision + recall)
    return( recall, precision, f1_score )



def main():
    parser = argparse.ArgumentParser(
        description="This script compute precision and recall of predictions based on taxonomic nodes comparison"
    )
    parser.add_argument(
        "--dataset",
        type=str,
        help="Dataset complexity type: high or low",
        required=True,
    )
    parser.add_argument(
        "--synth_data",
        type=str,
        help="Synthetic community dataset synth_[1..20]",
        required=True,
    )
    args = parser.parse_args()

    prefix_path = '/hps/nobackup/rdf/metagenomics/service-team/projects/finding-pheno/shallow_shotgun/'

    ## Catalogue lineages file
    catalogue_metadata = prefix_path + 'metadata_files/genomes-all_metadata.tsv'

    ### Files paths to results
    expected_covs = prefix_path + 'synthetic_communities/taxo_test/' + args.dataset + '_comp/coverage_tables/concat_' + args.dataset +'_cov.txt'

    bwa_mock_path = prefix_path + 'alignments/' + args.dataset + '/' + args.synth_data
    bwa_reps_all = bwa_mock_path + '/' + args.dataset + '_all.tsv'
    bwa_reps_best = bwa_mock_path + '/' + args.dataset + '_best.tsv'
    bwa_reps_unique = bwa_mock_path + '/' + args.dataset + '_unique.tsv'

    sm_mock_path = prefix_path + 'sourmash_results/' + args.dataset + '/' +args.synth_data
    sm_reps_11 = sm_mock_path + '/' + args.dataset + '_' + args.synth_data + '_k11_species_uw.tsv'
    sm_reps_21 = sm_mock_path + '/' + args.dataset + '_' + args.synth_data + '_k21_species_uw.tsv'
    sm_reps_31 = sm_mock_path + '/' + args.dataset + '_' + args.synth_data + '_k31_species_uw.tsv'
    sm_reps_51 = sm_mock_path + '/' + args.dataset + '_' + args.synth_data + '_k51_species_uw.tsv'


    ### Calling functions
    ( catalogue_lineages ) = metadata_parser( catalogue_metadata )
    ( mock_nodes ) = mock_parser( args.synth_data, expected_covs )

    bwa_nodes = {}
    ( bwa_nodes ) = bwa_parser( bwa_reps_all, 'all', catalogue_lineages, bwa_nodes )
    ( bwa_nodes ) = bwa_parser( bwa_reps_best, 'best', catalogue_lineages, bwa_nodes )
    ( bwa_nodes ) = bwa_parser( bwa_reps_unique, 'unique', catalogue_lineages, bwa_nodes )

    sm_nodes = {}
    ( sm_nodes ) = sm_parser( sm_reps_11, 'kmer_11', sm_nodes )
    ( sm_nodes ) = sm_parser( sm_reps_21, 'kmer_21', sm_nodes )
    ( sm_nodes ) = sm_parser( sm_reps_31, 'kmer_31', sm_nodes )
    ( sm_nodes ) = sm_parser( sm_reps_51, 'kmer_51', sm_nodes )

    bwa_metrics_writer( args.dataset, args.synth_data, mock_nodes, bwa_nodes )
    sm_metrics_writer( args.dataset, args.synth_data, mock_nodes, sm_nodes )


if __name__ == "__main__":
    main()

