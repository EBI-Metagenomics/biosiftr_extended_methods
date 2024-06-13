#!/usr/bin/env python

import argparse
import os.path
import sys
import gzip

##### This script transforms SHOGUN taxonomy results from NCBI taxonomy to GTDB
##### Alejandra Escobar, EMBL-EBI
##### June 12, 2024



def metadata_parser( metadata_r207 ):
    ncbi_names_dict = {}
    ncbi_organism = {}
    silva_names_dict = {}
    silva_23s_dict = {}
    for meta_file in metadata_r207:
        with gzip.open(meta_file, "rt") as input_file:
            next( input_file )
            for line in input_file:
                l_line = line.rstrip().split('\t')
                gtdb_taxonomy = l_line[16]
                ncbi_org = l_line[62].replace(' ','_')
                ncbi_taxonomy = l_line[78]
                ncbi_species = ncbi_taxonomy.split(';')[-1].replace(' ','_')
                if len(ncbi_species) > 3:
                    if ncbi_species in ncbi_names_dict:
                        if gtdb_taxonomy not in ncbi_names_dict[ncbi_species]:
                            ncbi_names_dict[ncbi_species].append(gtdb_taxonomy)
                    else:
                        ncbi_names_dict[ncbi_species] = [gtdb_taxonomy]

                if ncbi_org in ncbi_organism:
                    if gtdb_taxonomy not in ncbi_organism[ncbi_org]:
                        ncbi_organism[ncbi_org].append(gtdb_taxonomy)
                else:
                    ncbi_organism[ncbi_org] = [gtdb_taxonomy]

                silva_org = l_line[105].split(';')[-1].replace(' ','_')
                if silva_org in silva_names_dict:
                    if gtdb_taxonomy not in silva_names_dict[silva_org]:
                        silva_names_dict[silva_org].append(gtdb_taxonomy)
                else:
                    silva_names_dict[silva_org] = [gtdb_taxonomy]
                
                silva_23s_org = l_line[37].split(';')[-1].replace(' ','_')
                if silva_23s_org in silva_23s_dict:
                    if gtdb_taxonomy not in silva_23s_dict[silva_23s_org]:
                        silva_23s_dict[silva_23s_org].append(gtdb_taxonomy)
                else:
                    silva_23s_dict[silva_23s_org] = [gtdb_taxonomy]

    return( ncbi_names_dict, ncbi_organism, silva_names_dict, silva_23s_dict )


def lca_finder( lineage_list ):
    if len(lineage_list) > 1:
        split_lineages = [lineage.split(';') for lineage in lineage_list]
        common_lineage = split_lineages[0]
        for lineage in split_lineages[1:]:
            for i in range(len(common_lineage)):
                if i >= len(lineage) or common_lineage[i] != lineage[i]:
                    common_lineage = common_lineage[:i]
                    break
        gtdb_name = ';'.join(common_lineage)
    else:
        gtdb_name = lineage_list[0]
    gtdb_name = gtdb_name.replace(' ','_')

    return( gtdb_name )


def shogun_parser( shogun_tax, ncbi_names_dict, ncbi_organism, silva_names_dict, silva_23s_dict, output_name ):
    accepted_domains = [
        "k__Bacteria",
        "k__Archaea"
    ]
    curated = {
        "s__Parageobacillus_thermoglucosidans":"s__Parageobacillus_thermoglucosidasius",
        "s__Atopobium_parvulum":"Lancefieldella_parvula",
	"s__Clostridiales_bacterium_oral_taxon_876":"s__Clostridiales_bacterium_oral_taxon_876_str._F0540",
	"s__[Clostridium]_ultunense":"s__Schnuerera_ultunensis",
	"s__Candidatus_Arthromitus_sp._SFB-mouse":"s__Candidatus_Arthromitus_sp._SFB-mouse-NL",
	"s__cyanobacterium_endosymbiont_of_Epithemia_turgida":"s__cyanobacterium_endosymbiont_of_Epithemia_turgida_isolate_EtSB_Lake_Yunoko",
	"s__[Clostridium]_asparagiforme":"s__Enterocloster_asparagiformis",
	"s__[Clostridium]_clariflavum":"s__Acetivibrio_clariflavus",
	"s__Loktanella_cinnabarina":"s__Limimaricola_cinnabarinus",
	"s__Aster_yellows_witches'-broom_phytoplasma":"s__Aster_yellows_witches'-broom_phytoplasma_AYWB",
	"s__Anaerococcus_provenciensis":"s__Anaerococcus_provencensis",
	"s__Megasphaera_genomosp._type_1":"s__Megasphaera_genomosp._type_1_str._28L",
	"s__Propionibacterium_sp._oral_taxon_192":"s__Propionibacterium_sp._oral_taxon_192_str._F0372",
	"s__[Clostridium]_purinilyticum":"s__Gottschalkia_purinilytica",
	"s__Streptomyces_sp._S10(2016)":"s__Streptomyces_qaidamensis",
	"s__Bradyrhizobium_genosp._SA-4":"s__Bradyrhizobium_genosp._SA-4_str._CB756",
	"s__Arthrobacter_sp._IHBB_11108":"s__Psychromicrobium_lacuslunae",
	"s__Oribacterium_sp._oral_taxon_078":"s__Oribacterium_sp._oral_taxon_078_str._F0263",
	"s__Cyanothece_sp._PCC_8801":"s__Rippkaea_orientalis_PCC_8801",
	"s__Actinomyces_odontolyticus":"s__Actinomyces_odontolyticus_F0309",
	"s__Olsenella_sp._oral_taxon_809":"s__Olsenella_sp._oral_taxon_809_str._F0356",
	"s__Coprobacillus_sp._D6":"s__Coprobacillus_sp._29_1",
	"s__Anaerobacillus_macyae":"s__Alkalihalobacillus_macyae",
	"s__Desulfotomaculum_alkaliphilum":"s__Desulfohalotomaculum_alkaliphilum_DSM_12257",
	"s__[Clostridium]_aerotolerans":"s__Desulfovibrio_aerotolerans",
	"s__Lactobacillus_thailandensis":"s__Lactobacillus_thailandensis_DSM_22698_=_JCM_13996",
	"s__Lactobacillus_alimentarius":"s__Lactobacillus_alimentarius_DSM_20249",
	"s__Desulfotomaculum_acetoxidans":"s__Desulfofarcimen_acetoxidans",
	"s__Lactobacillus_vaginalis":"s__Lactobacillus_vaginalis_DSM_5837_=_ATCC_49540",
	"s__Methylobacterium_populi":"s__Methylorubrum_populi",
	"s__Ruminococcus_faecis":"s__Mediterraneibacter_faecis",
	"s__Bacillus_thermotolerans":"s__Bacillaceae_bacterium_MTCC_8252",
	"s__Lactobacillus_antri":"s__Lactobacillus_antri_DSM_16041",
	"s__Bifidobacterium_gallinarum":"s__Bifidobacterium_pullorum_subsp._gallinarum",
	"s__Lactobacillus_uvarum":"s__Lactobacillus_uvarum_DSM_19971",
	"s__Lactobacillus_hordei":"s__Lactobacillus_hordei_DSM_19519",
	"s__Chryseobacterium_caeni":"s__Chryseobacterium_caeni_DSM_17710",
	"s__Lactobacillus_graminis":"s__Lactobacillus_graminis_DSM_20719",
	"s__Bacillus_bogoriensis":"s__Bacillus_bogoriensis_ATCC_BAA-922",
	"s__Actinomyces_europaeus":"s__Actinomyces_europaeus_ACS-120-V-Col10b",
	"s__Lysinibacillus_manganicus":"s__Lysinibacillus_manganicus_DSM_26584",
	"s__Enterorhabdus_caecimuris":"s__Enterorhabdus_caecimuris_B7",
	"s__Atopobium_vaginae":"s__Fannyhessea_vaginae",
	"s__Escherichia_vulneris":"s__Pseudescherichia_vulneris",
	"s__Actinomyces_turicensis":"s__Actinomyces_turicensis_ACS-279-V-Col4",
	"s__Mycobacterium_aromaticivorans":"s__Mycobacterium_aromaticivorans_JS19b1_=_JCM_16368",
	"s__Bacillus_solani":"s__Cytobacillus_solani",
	"s__Tetrasphaera_elongata":"s__Tetrasphaera_elongata_Lp2",
	"s__Desulfovibrio_alkalitolerans":"s__Desulfovibrio_alkalitolerans_DSM_16529",
	"s__Pedobacter_oryzae":"s__Pedobacter_oryzae_DSM_19973",
	"s__Lysinimicrobium_mangrovi":"s__Demequina_mangrovi",
	"s__[Clostridium]_termitidis":"s__Ruminiclostridium_cellobioparum_subsp._termitidis",

    }
    hardcoded = {
        "s__Bacillus_alveayuensis":"d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Aeribacillaceae;g__Bacillus_CB;s__Bacillus_CB_alveayuensis",
        "s__Rheinheimera_texasensis":"d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Alteromonadaceae;g__Pararheinheimera;s__Pararheinheimera_texasensis"
    }
    aggregated = {}
    singleton_counter = 0
    with open(shogun_tax, 'r') as input_file, open(output_name, 'w') as output_file:
        header = input_file.readline().strip()
        for line in input_file:
            line_l = line.rstrip().split('\t')
            sho_name = line_l.pop(0)

            values = [float(item) for item in line_l]
            val_sum = sum(values)
            if val_sum > 2:
                domain = sho_name.split(';')[0]
                if domain in accepted_domains:
                    current_values = [float(item) for item in line_l]
                    species = sho_name.split(';')[-1]

                    if species in curated:
                        species = curated[species]

                    if species in ncbi_names_dict:
                        gtdb_name = lca_finder( ncbi_names_dict[species] )
                    elif species.replace('s__','') in ncbi_organism:
                        gtdb_name = lca_finder( ncbi_organism[species.replace('s__','')] )
                    elif species.replace('s__','') in silva_names_dict:
                        gtdb_name = lca_finder( silva_names_dict[species.replace('s__','')] )
                    elif species.replace('s__','') in silva_23s_dict:
                        gtdb_name = lca_finder( silva_23s_dict[species.replace('s__','')] )
                    elif species in hardcoded:
                        gtdb_name = hardcoded[species]
                    else:
                        strain = species.split('_')[-1]
                        genomes_list = []
                        for ncbi_name in ncbi_names_dict:
                            ncbi_strain = ncbi_name.split('_')[-1]
                            if strain == ncbi_strain:
                                genomes_list = genomes_list + ncbi_names_dict[ncbi_name]
                        if len(genomes_list) > 0:
                            gtdb_name = lca_finder( genomes_list )
                        else:
                            print("No GTDB annotation for NCBI species: "+species)

                    if gtdb_name in aggregated:
                        saved_values = aggregated[gtdb_name]
                        added_values = [x + y for x, y in zip(saved_values, current_values)]
                        aggregated[gtdb_name] = added_values
                    else:
                        aggregated[gtdb_name] = current_values

            else:
                singleton_counter+=1

    print('Discarded singletons and doubletons: '+str(singleton_counter))

    # Writing the output
    with open(output_name, 'w') as output_file:
        output_file.write( header + '\n' )
        for lineage in aggregated:
            values = [str(item) for item in aggregated[lineage]]
            output_file.write( lineage + '\t' + '\t'.join(values)+ '\n')


def main():
    parser = argparse.ArgumentParser(
            description="This script transforms SHOGUN taxonomy results from NCBI names to GTDB r207"
    )
    parser.add_argument(
        "--shogun_tax",
        type=str,
        help="Integrated table of taxonomy generated by SHOGUN",
        required=True,
    )
    parser.add_argument(
        "--metadata_r207",
        nargs=2,
        help="List of GTDB input files: ar53_metadata_r207.tsv.gz and bac120_metadata_r207.tsv.gz",
        required=True,
    )
    parser.add_argument(
       "--output",
        type=str,
        help="Name of the output file",
        required=True,
    )
    args = parser.parse_args()
    
    ### Calling functions
    ( ncbi_names_dict, ncbi_organism, silva_names_dict, silva_23s_dict ) = metadata_parser( args.metadata_r207 )
    shogun_parser( args.shogun_tax, ncbi_names_dict, ncbi_organism, silva_names_dict, silva_23s_dict, args.output )


if __name__ == "__main__":
    main()

