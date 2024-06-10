#!/usr/bin/env python

import argparse
import os.path
import pysam
import sys
from Bio import SeqIO

##### This script generates a table of the expected and observed relative abundance for the shallow-shotgun paper
##### Alejandra Escobar, EMBL-EBI
##### Dec 19, 2023


def exp_parser( ref_abun, synth_prefix ):
    exp_dict = {}
    with open(ref_abun, 'r') as input_file:
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
            if dataset == synth_prefix:
                genome = genome.replace('.fa','')
                #mock_lineage = mock_lineage.replace('p__Bacillota','p__Firmicutes')
                exp_dict[mock_lineage] = assem_abundance
    return( exp_dict )


def obs_parser( obs_table, exp_dict, synth_prefix ):
    obs_abun = {}
    with open(obs_table, 'r') as input_file:
        next(input_file)
        for line in input_file:
            lineage_g,rel_abun = line.rstrip().split('\t')
            lineage = lineage_g.replace('g__;','').replace('s__;','').split(';')
            lineage.pop(-1)
            lineage = ';'.join(lineage)

            if lineage in obs_abun:
                sum_abun = float(obs_abun[lineage]) + float(rel_abun)
                obs_abun[lineage] = sum_abun
            else:
                obs_abun[lineage] = float(rel_abun)

    
    exp_used = []
    with open('relab_obs_exp.tsv', 'w') as output_file:
        output_file.write("\t".join([
            'dataset',
            'species',
            'exp_relab',
            'obs_relab'])+'\n')
        for lineage in obs_abun:
            obs_relab = obs_abun[lineage]
            if lineage in exp_dict:
                exp_used.append(lineage)
                exp_relab = exp_dict[lineage]
            else:
                exp_relab = '0'
            output_file.write("\t".join([
                synth_prefix,
                lineage,
                exp_relab,
                str(obs_relab)])+'\n')

        for lineage in exp_used:
            del exp_dict[lineage]

        for lineage in exp_dict:
            exp_abun = exp_dict[lineage]
            output_file.write("\t".join([
                synth_prefix,
                lineage,
                exp_abun,
                '0'])+'\n')


def main():
    parser = argparse.ArgumentParser(
        description="This script generates a table of the expected and observed relative abundance for the shallow-shotgun paper"
    )
    parser.add_argument(
        "--sp_relab",
        type=str,
        help="Relative abundance table at species level (species_relab.tsv)",
        required=True,
    )
    parser.add_argument(
        "--dataset_type",
        type=str,
        help="Either low or high",
        required=True,
    )
    parser.add_argument(
        "--synth_prefix",
        type=str,
        help="Synhtetic dataset synth_[1..20]",
        required=True,
    )
    args = parser.parse_args()


    ref_abun = '/hps/nobackup/rdf/metagenomics/service-team/projects/finding-pheno/shallow_shotgun/synthetic_communities/taxo_test/' + args.dataset_type + '_comp/coverage_tables/concat_' + args.dataset_type + '_cov.txt'

    ### Calling functions
    ( exp_dict ) = exp_parser( ref_abun, args.synth_prefix )
    obs_parser( args.sp_relab, exp_dict, args.synth_prefix )


if __name__ == "__main__":
    main()

