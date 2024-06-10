#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob
import gzip


##### This script calculates the coverage per genome of a microbial community from the metagenomic raw-reads and the assemblies used for InSilicoSeq input
##### Alejandra Escobar, EMBL-EBI
##### Jul 6, 2023


def fq_parser( fq_in ):
    reads_counts = {}
    with gzip.open(fq_in, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            contig_id = str(record.id).split('_')
            contig_id.pop(-1)
            contig_id.pop(-1)
            contig_id = '_'.join(contig_id)
            if contig_id in reads_counts:
                reads_counts[contig_id]+=1
            else:
                reads_counts[contig_id]=1

    return( reads_counts )


def abun_parser( abun_in ):
    rel_ab = {}
    head, tail = os.path.split(abun_in)
    dataset = tail.replace('_abundance.txt','') 
    with open(abun_in, 'r') as file_in:
        for line in file_in:
            genome,abundance = line.rstrip().split("\t")
            dir_name, genome = genome.split('/')
            genome = genome.replace('.fasta','')
            rel_ab[genome] = abundance

    return( dataset, rel_ab )


def fasta_parser( assemblies ):
    genomes_size = {}
    contig_ids = {}
    for assembly in assemblies:
        head, tail = os.path.split(assembly)
        #tail = tail.replace('.fa','')
        genome_len = 0
        contig_ids[tail] = []
        for record in SeqIO.parse(assembly, "fasta"):
            contig_id = str(record.id)
            genome_len = genome_len + int(len(str(record.seq)))
            contig_ids[tail].append(contig_id)
        genomes_size[tail] = genome_len

    return( genomes_size, contig_ids )

   
def table_writer( reads_counts, rel_ab, genomes_size, contig_ids, dataset ):
    with open(dataset+'_genomes_cov.txt', 'w') as file_out:
        file_out.write('\t'.join([
            'genome',
            'dataset',
            'total_reads',
            'genome_size',
            'assem_abundance',
            'genome_coverage'
        ]) + '\n')
        for assembly in rel_ab:
            assem_abundance = rel_ab[assembly]
            genome_len = genomes_size[assembly]
            total_reads = 0
            for contig in contig_ids[assembly]:
                if contig in reads_counts:
                    total_reads = total_reads+reads_counts[contig]
            coverage = ( float(300) * float(total_reads) ) / float(genome_len)
            coverage = round(coverage, 3)
            file_out.write('\t'.join([
                assembly, 
                dataset, 
                str(total_reads), 
                str(genome_len),
                assem_abundance, 
                str(coverage)
            ]) + '\n')

def main():
    parser = argparse.ArgumentParser(
        description="This script calculates the coverage per genome of a microbial community from the metagenomic raw-reads and the assemblies used for InSilicoSeq input. Please provide the relevant input files"
    )
    parser.add_argument(
        "--read_1",
        type=str,
        help="Raw reads output by InSilicoSeq, pair one only (file_1.fq.gz)",
        required=True,
    )
    parser.add_argument(
        "--assem_list",
        nargs="*", 
        help="List of assemblies used as input for InSilicoSeq", 
        required=True,
    )
    parser.add_argument(
        "--abun_table",
        type=str,
        help="Abundance profile output by InSilicoSeq",
        required=True,
    )
    args = parser.parse_args()

    ## Calling functions
    # Counting raw-reads
    reads_counts = fq_parser( args.read_1 )

    # Saving the relative abundance values
    ( dataset, rel_ab ) = abun_parser( args.abun_table )

    # Saving the genomes size and contig ids
    ( genomes_size, contig_ids ) = fasta_parser( args.assem_list )
    
    # Generating the coverage table
    table_writer( reads_counts, rel_ab, genomes_size, contig_ids, dataset )


if __name__ == "__main__":
    main()


