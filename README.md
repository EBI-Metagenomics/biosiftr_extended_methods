# shallow_shotgun_paper
This repo contains the scripts and tables generated to optimise and validate the [MGnify shallowmapping tool](https://github.com/EBI-Metagenomics/shallowmapping). Described in the following publication:


// TODO Add publication reference

Leveraging MGnify Genomic Catalogues for Inferring Metabolic Potential in Shallow-Shotgun Sequenced Samples. Alejandra Escobar-Zepeda, Matti Ruuskanen, Martin Beracochea, Lorna Richardson, Robert D. Finn, Leo Lahti.

<p align="center" width="100%">
   <img src="visuals/abstract.png" width="100%"/>
</p>

## Contents
- [Section 1 . Shallow-mapping tool optimisation](#sec1)
  1. Synthetic communities design
  2. Taxonomic profile prediction power
  3. Functional annotation benchmark
- [Section 2. Pipeline validation on real data](#sec2)
  1. Comparative metagenomics of KOs presence/absence


<a name="sec1"></a>
## Section 1 . Shallow-mapping tool optimisation
### 1. Synthetic communities design

To optimise the parameters for [bwamem2](https://github.com/bwa-mem2/bwa-mem2) and [Sourmash](https://github.com/sourmash-bio/sourmash) mapping tools performance in the [MGnify Shallow-mapping tool](https://github.com/EBI-Metagenomics/shallowmapping), we generated synthetic microbial communities according to the following schema.

<p align="center" width="100%">
   <img src="visuals/synthetic_shallow.png" width="100%"/>
</p>

[InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) was used to generate each of the synthetic communities from [MetaChic](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FHPJH5) chicken-gut genomes catalogue using representative genomes only. The genome contigs were renamed to make it easier to track back the origin of the synthetic reads after mapping

```bash
# Renaming the contigs of all representatives
assembly_rename.py rep_genome.fa rep_genome

# Preparing genomes for tools optimisation (genomes with annotation at species level present in the MGnify catalogue)
cd taxonomy/rich_50
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../known_realpath_genomes.txt | shuf -n50); do (ln -s $genome .); done); done

cd taxonomy/rich_500
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../known_realpath_genomes.txt | shuf -n500); do (ln -s $genome .); done); done

# Preparing genomes for functional benchmarking (all genomes in MetaChick catalogue)
cd function/rich_50
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../all_realpath_genomes.txt | shuf -n50); do (ln -s $genome .); done); done

cd function/rich_500
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../all_realpath_genomes.txt | shuf -n500); do (ln -s $genome .); done); done

# Generating shallow-shotgun raw-reads
for num in {1..20}; do (iss generate --draft synth_$num/*.fa --model novaseq --output raw_reads/synth_$num --cpus 4 --n_reads 2M --compress

```

The ground truth to benchmark functional prediction power was generated using [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) directly from MetaChick MAGs used to build the synthetic communities.

```bash
# Generating functional annotation for each MetaChick MAG
emapper.py -i rep_genome.fa --itype genome --database eggnog/data/eggnog.db --dmnd_db eggnog/data/eggnog_proteins.dmnd --data_dir eggnog/data/ -m diamond --no_file_comments --cpu 16 --dbmem

# Integrating eggnog functions of synthetic communities (emapper evalue threshold = 1E-10)
synth_functions_integrator.py --eggnog metachick_representative_mags/emapper_results --dataset_type low
```

### 2. Taxonomic profile prediction power

We used the MGnify chicken-gut catalogue representative genomes as a reference for mapping. The tools were run with the corresponding parameters and the outputs were processed as described below.

```bash
## bwamem2
# Indexing MGnify representative genomes for bwamem2
bwa-mem2 index -p bwa_reps.fa concat_mgnify_reps.fa

# Generating bwamem2 alignments
for num in {1..20}; do (cd synth_$num && bwa-mem2 mem -M -t 8 bwa_reps.fa synth_$num\_R1.fastq.gz synth_$num\_R2.fastq.gz | samtools view -@ 8 -F256 -F4 -uS - | samtools sort -@ 8 -O bam - -o sort_filt_reps.bam && samtools index sort_filt_reps.bam); done

# QC (ani and cov) filtering and processing bwamem2 alignment types
bam2cov.py --bwa_bam sort_filt_reps.bam --prefix low

# A step to transform genomes relative abundance to species relative abundance. The genomes-all_metadata.tsv file is available at the chicken git catalogue ftp site: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/chicken-gut/v1.0.1/genomes-all_metadata.tsv
bwa_genome2species.py --genomes_relab u_relab_1.tsv --metadata genomes-all_metadata.tsv --output low_species_relab.tsv

# Expected vs bwamem2 observed relative abundance
bwa_exp_obs_relab.py --sp_relab low_species_relab.tsv --dataset_type low --synth_prefix synth_$num


## Sourmash
# Indexing the representative species genomes for Sourmash db kmer = 11, 21, 31, 51
sourmash sketch dna -p scaled=1000,k=$kmer --name "$genome" $genome.fna -o $genome.fna.sig
sourmash index chicken_gut_k${kmer}.sbt.zip *sig

# Sketching the raw-reads for Sourmash
sourmash sketch dna -p k=$kmer,abund synth_$num\_R1.fastq.gz synth_$num\_R2.fastq.gz --name synth_$num\_$kmer -o synth_$num\_$kmer.zip

# Running Sourmash gather on representative species dbs
sourmash gather -o low_synth_$num\_gather-k$kmer.csv synth_$num\_$kmer.zip reps_chicken_gut_k$kmer.sbt.zip

# Editing the Sourmash CSV output to add lineage and keep relevant info. The genomes-all_metadata.tsv file is available at the chicken git catalogue ftp site: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/chicken-gut/v1.0.1/genomes-all_metadata.tsv 
sm_genome2species.py --sm_csv low_synth_$num\_gather-k$kmer.csv --metadata genomes-all_metadata.tsv --output low_synth_$num\_k$kmer\_species_uw.tsv

# Expected vs Sourmash observed relative abundance
sm_exp_obs_relab.py --sp_relab low_synth_$num\_k51_species_uw.tsv --dataset_type low --synth_data synth_$num


## Processing bwamem2 and sourmash results
# Calculating expected coverage. Abundance tables are generated by InSilicoSeq
coverage_calc.py --read_1 raw_reads/synth_$num\_R1.fastq.gz --assem_list synth_$num/*.fa --abun_table synth_$num\_abundance.txt

# Calculating performance metrics (F1, precision, recall)
tree_precrecal.py --dataset low --synth_data synth_$num

```



