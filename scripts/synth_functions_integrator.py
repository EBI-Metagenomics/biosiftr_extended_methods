#!/usr/bin/env python

import argparse
import os.path
import gzip
import sys
from Bio import SeqIO

##### This script generates the functioal annotation tables for synthetic communities of the shallow-shotgun paper
##### Alejandra Escobar, EMBL-EBI
##### Jan 11, 2024


def pfam_parser( pfam_data ):
    pfam_desc = {}
    with gzip.open(pfam_data, 'rt', encoding='utf-8') as input_file:
        entry = []
        for line in input_file:
            line = line.strip()
            if line == '//' and entry:
                entry_dict = {}
                for entry_line in entry:
                    if entry_line.startswith('#=GF ID'):
                        entry_dict['ID'] = entry_line.split()[-1]
                    elif entry_line.startswith('#=GF AC'):
                        entry_dict['AC'] = entry_line.split()[-1].split('.')[0]
                if 'ID' in entry_dict and 'AC' in entry_dict:
                    pfam_desc[entry_dict['ID']] = entry_dict['AC']
                entry = []
            else:
                entry.append(line)
    return( pfam_desc )


def exp_parser( ref_abun ):
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
            genome = genome.replace('.fa','').replace('renamed_','')
            if dataset in exp_dict:
                exp_dict[dataset].append(genome)
            else:
                exp_dict[dataset] = [genome]
    return( exp_dict )


def eggnog_parser( exp_dict, eggnog_path, pfam_desc ):
    eggnog_dict = {}
    all_kos = []
    all_pfams = []
    for dataset in exp_dict:
        for genome in exp_dict[dataset]:
            eggnog_file = eggnog_path+'/'+genome+'_out.emapper.annotations'
            with open(eggnog_file, 'r') as input_file:
                next(input_file)
                for line in input_file:
                    l_line = line.rstrip().split('\t')

                    evalue = float(l_line[2])
                    if evalue <= 1e-10:
                        kegg_ko = l_line[11].replace('ko:','')
                        if kegg_ko != '-':
                            if ',' in kegg_ko:
                                ko_list = kegg_ko.split(',')
                            else:
                                ko_list = [kegg_ko]

                            for ko_id in ko_list:
                                all_kos.append(ko_id)
                                composite_key = ( dataset, ko_id )
                                if composite_key in eggnog_dict:
                                    eggnog_dict[composite_key] += 1
                                else:
                                    eggnog_dict[composite_key] = 1

                        pfam_key = (genome, 'pfam')
                        pfams_raw = l_line[20]
                        if pfams_raw != '-':
                            pfam_ids_list = []
                            if ',' in pfams_raw:
                                for pfam in pfams_raw.split(','):
                                    if pfam in pfam_desc:
                                        pfam_id = pfam_desc[pfam]
                                        pfam_ids_list.append(pfam_id)
                            else:
                                if pfams_raw in pfam_desc:
                                    pfam_ids_list.append(pfam_desc[pfams_raw])
                            if len(pfam_ids_list)> 0:
                                for pfam_annot in pfam_ids_list:
                                    all_pfams.append(pfam_annot)
                                    composite_key = ( dataset, pfam_annot )
                                    if composite_key in eggnog_dict:
                                        eggnog_dict[composite_key] += 1
                                    else:
                                        eggnog_dict[composite_key] = 1

    all_kos = list(set(all_kos))
    all_pfams = list(set(all_pfams))
    return( eggnog_dict, all_kos, all_pfams )


def table_writer( dataset_type, eggnog_dict, all_kos, all_pfams ):
    dataset_list = []
    for num in range(1,21):
        value = 'synth_'+str(num)
        dataset_list.append(value)

    # Generating the kegg counts table
    with open(dataset_type+'_kegg_counts.tsv', 'w') as kegg_out:
        header = ['ko_id'] + dataset_list
        kegg_out.write('\t'.join(header)+'\n')
        for ko_id in all_kos:
            to_print = []
            to_print.append(ko_id)
            for dataset in dataset_list:
                composite_key = (dataset, ko_id)
                if composite_key in eggnog_dict:
                    to_print.append(str(eggnog_dict[composite_key]))
                else:
                    to_print.append('0')
            kegg_out.write('\t'.join(to_print)+'\n')

    # Generating the Pfam counts table
    with open(dataset_type+'_pfam_counts.tsv', 'w') as pfam_out:
        header = ['pfam_id'] + dataset_list
        pfam_out.write('\t'.join(header)+'\n')
        for pfam_id in all_pfams:
            to_print = []
            to_print.append(pfam_id)
            for dataset in dataset_list:
                composite_key = (dataset, pfam_id)
                if composite_key in eggnog_dict:
                    to_print.append(str(eggnog_dict[composite_key]))
                else:
                    to_print.append('0')
            pfam_out.write('\t'.join(to_print)+'\n')




def main():
    parser = argparse.ArgumentParser(
        description="This script generates the functioal annotation tables for synthetic communities of the shallow-shotgun paper"
    )
    parser.add_argument(
        "--eggnog",
        type=str,
        help="EggNOG results path",
        required=True,
    )
    parser.add_argument(
        "--dataset_type",
        type=str,
        help="Dataset complexity type: low or high",
        required=True,
    )
    args = parser.parse_args()

    ref_abun = '/hps/nobackup/rdf/metagenomics/service-team/projects/finding-pheno/shallow_shotgun/synthetic_communities/func_test/' + args.dataset_type + '_comp/coverage_tables/concat_' + args.dataset_type + '_cov.txt'

    ### Calling functions
    ( pfam_desc ) = pfam_parser( '/hps/nobackup/rdf/metagenomics/service-team/ref-dbs/dram/Pfam-A.hmm.dat.gz' )
    ( exp_dict ) = exp_parser( ref_abun )
    ( eggnog_dict, all_kos, all_pfams ) = eggnog_parser( exp_dict, args.eggnog, pfam_desc )
    table_writer( args.dataset_type, eggnog_dict, all_kos, all_pfams )


if __name__ == "__main__":
    main()

