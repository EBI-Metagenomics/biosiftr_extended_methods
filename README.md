# BioSIFTR methods

This repo contains the scripts and tables generated to optimise and validate the [MGnify BioSIFTR tool](https://github.com/EBI-Metagenomics/biosiftr) described in the following publication:

Biome-specific genome catalogues reveal functional potential of shallow sequencing. Alejandra Escobar-Zepeda, Matti O. Ruuskanen, Martin Beracochea, Jennifer Lu, Dattatray Mongad, Lorna Richardson, Robert D. Finn, Leo Lahti

<p align="center" width="100%">
   <img src="visuals/visual_abstract_v2.png" width="100%"/>
</p>

## Contents
- [Section 1. The BioSIFTR tool optimisation](#sec1)
  1. Synthetic communities design
  2. Taxonomic profile prediction power
  3. Generation of pangenome tables
  4. Functional annotation benchmark
- [Section 2. Pipeline validation on real data](#sec2)
  1. Functional annotation and comparative metagenomics analysis
  2. Taxonomic annotation


<a name="sec1"></a>
## Section 1. The BioSIFTR tool optimisation
### 1. Synthetic communities design

To optimise the parameters for [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) and [Sourmash](https://github.com/sourmash-bio/sourmash) mapping tools' performance in the [MGnify BioSIFTR tool](https://github.com/EBI-Metagenomics/biosiftr), we generated synthetic microbial communities according to the following schema.

<p align="center" width="100%">
   <img src="visuals/synthetic_shallow.png" width="100%"/>
</p>

[InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) was used to generate each of the synthetic communities from [MetaChic](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FHPJH5) chicken-gut genomes catalogue using representative genomes only. The genome contigs were renamed to make it easier to trace back the origin of the synthetic reads after mapping. For the human-gut, we generated an external catalogue from a mix of isolates, MAGs and SAGs. Accessions are available in the supplementary material of the actual publication.

```bash
# Renaming the contigs of all representatives
assembly_rename.py rep_genome.fa rep_genome

# Preparing genomes for tools optimisation (genomes with annotation at species level present in the MGnify catalogue)
cd taxonomy/rich_50
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../known_realpath_genomes.txt | shuf -n50); do (ln -s $genome .); done); done

cd taxonomy/rich_500
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../known_realpath_genomes.txt | shuf -n500); do (ln -s $genome .); done); done

# Preparing genomes for functional benchmarking (all genomes in the external catalogue)
cd function/rich_50
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../all_realpath_genomes.txt | shuf -n50); do (ln -s $genome .); done); done

cd function/rich_500
for num in {1..20}; do (cd synth_$num && for genome in $(cat ../all_realpath_genomes.txt | shuf -n500); do (ln -s $genome .); done); done

# Generating shallow-shotgun raw reads
for num in {1..20}; do (iss generate --draft synth_$num/*.fa --model novaseq --output raw_reads/synth_$num --cpus 4 --n_reads 2M --compress); done

```

The ground truth to benchmark the functional prediction power was generated using [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) directly from the external catalogue of MAGs used to build the synthetic communities.

```bash
# Generating functional annotation for each MetaChick MAG
emapper.py -i rep_genome.fa --itype genome --database eggnog/data/eggnog.db --dmnd_db eggnog/data/eggnog_proteins.dmnd --data_dir eggnog/data/ -m diamond --no_file_comments --cpu 16 --dbmem

# Integrating eggnog functions of synthetic communities (emapper e-value threshold = 1E-10)
synth_functions_integrator.py --eggnog metachick_representative_mags/emapper_results --dataset_type low
```

### 2. Taxonomic profile prediction power
We used the MGnify chicken-gut catalogue representative genomes as a reference database. Mapping reads generated with bwamem2 were classified as follows:
- unique: reads mapping to one reference only (no XA tag on the bam file)
- all: unique reads plus all extra matches (every reference on the XA:Z tag having a minimum of 60% of the read length of exact matches is counted in)
- best: unique reads plus reads with XA:Z tag. Only the match reported in the third column of the BAM file is considered.

```bash
## bwamem2
# Indexing MGnify representative genomes for bwamem2
bwa-mem2 index -p bwa_reps.fa concat_mgnify_reps.fa

# Generating bwamem2 alignments
for num in {1..20}; do (cd synth_$num && bwa-mem2 mem -M -t 8 bwa_reps.fa synth_$num\_R1.fastq.gz synth_$num\_R2.fastq.gz | samtools view -@ 8 -F256 -F4 -uS - | samtools sort -@ 8 -O bam - -o sort_filt_reps.bam && samtools index sort_filt_reps.bam); done

# QC (ani and cov) filtering and processing bwa-mem2 alignment types. This script generates results for all, best, and unique mapping reads.
bam2cov.py --bwa_bam sort_filt_reps.bam --prefix low

# A step to transform genomes' relative abundance to species' relative abundance. The genomes-all_metadata.tsv file is available at the chicken git catalogue ftp site: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/chicken-gut/v1.0.1/genomes-all_metadata.tsv
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


## Processing bwamem2 and Sourmash results
# Calculating expected coverage. Abundance tables are generated by InSilicoSeq
coverage_calc.py --read_1 raw_reads/synth_$num\_R1.fastq.gz --assem_list synth_$num/*.fa --abun_table synth_$num\_abundance.txt

# Calculating performance metrics (F1, precision, recall)
tree_precrecal.py --dataset low --synth_data synth_$num

```

Performance metrics tables generated in this section are available in the [data/optimisation
/taxonomy/](https://github.com/EBI-Metagenomics/shallow_shotgun_paper/tree/main/data/optimisation/taxonomy) directory of this repo. 


### 3. Generation of pangenome tables

The optimisation of taxonomic annotation allows the accurate detection of species clusters in the MGnify genomes catalogue. Functional inference in the BioSIFTR pipeline is made through the annotation transference of the pangenome (or core genes). Tables of functions are generated from the tables available in the FTP site of the [MGnify genome catalogues](https://www.ebi.ac.uk/metagenomics/browse/genomes). This can be done by downloading the catalogue and then processing the files locally. The commands below were used in the local copy of the catalogue at the EBI cluster Codon. KEGG modules completeness tables were computed using the [kegg-pathways-completeness-tool](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool).

``` bash
# Generate the fasta files of genes lacking functional annotation
pangenomeDB_builder.py \
  --metadata genomes-all_metadata.tsv \
  --mode pre \
  --pfam_dat Pfam-A.hmm.dat.gz
 
# Launching the eggNOG annotation tool
for genome in $(ls *.fasta | sed 's/_accessory.fasta//'); do (emapper.py -i $genome\_accessory.fasta --itype CDS --translate --database eggnog/data/eggnog.db --dmnd_db data/eggnog_proteins.dmnd --data_dir eggnog/data/ -m diamond --no_file_comments --cpu 16 --dbmem -o $genome\_out); done
  
# Once all annotations have been completed successfully, move the annotation results into the emapper results directory and remove the intermediate files.
mkdir emapper_results
mv *_out.emapper.annotations emapper_results1
rm *.fasta *.hits *.seed_orthologs
 
# Parsing the annotation files and generating pre-computed profiles at the pangenome level
pangenomeDB_builder.py \
  --metadata genomes-all_metadata.tsv \
  --mode post \
  --pfam_dat Pfam-A.hmm.dat.gz
 
# Move the `*.tsv` files to a directory called `functional_profiles_DB`
mkdir functional_profiles_DB && mv *_clstr.tsv functional_profiles_DB
rm -r emapper_results

# Create the KEGG pathways database directory
mkdir kegg_completeness_DB
 
# Generating a KOs list from functional annotation tables for core and pangenome
for genome in $(ls functional_profiles_DB | cut -d'_' -f1); do (echo "processing core: " $genome && cat functional_profiles_DB/$genome\_clstr.tsv | sed '/^#/d' | awk -F '\t' '{if ($9 ~ /true/) print $6 }' | sed '/^-$/d' | tr "," "\n" | sort -u | tr "\n" "," | sed 's/,$/\n/' > kegg_completeness_DB/$genome\_core_kos.txt ); done

for genome in $(ls functional_profiles_DB | cut -d'_' -f1); do (echo "processing pan: " $genome && cat functional_profiles_DB/$genome\_clstr.tsv | sed '/^#/d' | cut -f6 | sed '/^-$/d' | tr "," "\n" | sort -u | tr "\n" "," | sed 's/,$/\n/' > kegg_completeness_DB/$genome\_pan_kos.txt ); done
 
# Running the KEGG pathways completeness tool
cd kegg_completeness_DB
for input in $(ls | cut -d'_' -f1,2 ); do (run_pathways.sh -l $input\_kos.txt -o $input); done
 
# Integrating the outputs
for genome in $(ls *.tsv | cut -d'_' -f1 | sort -u ); do (echo "processing: "$genome && keggcomp_DB.py --core $genome\_core.summary.kegg_pathways.tsv --pan $genome\_pan.summary.kegg_pathways.tsv --output $genome\_clstr); done
rm *.log *.txt *.kegg_pathways.tsv *.kegg_contigs.tsv

```

Pangenome tables generated by MGnify team are available in the FTP site of the corresponding MGnify genomes catalogue.

#### 3.1. Processing custom genome catalogues
Custom MAG catalogues generated using a minimal version of the [MGnify genomes-catalogue-pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) can be used to build bioSIFTR databases. The starting point is to have in the same directory (`genomes_dir`) the fasta files of the genomes to be included in the catalogue. You will need to install the following tools locally. We recommend using the publicly available Docker images:

- [CheckM v1.2.2](https://github.com/Ecogenomics/CheckM)
- [dRep v3.4.2](https://github.com/MrOlm/drep)
- [GTDBtk v2.4.0](https://github.com/Ecogenomics/GTDBTk)
- [PROKKA 1.14.6](https://github.com/tseemann/prokka)
- [EggNOG mapper v2.1.11](https://github.com/eggnogdb/eggnog-mapper)
- [Panaroo v1.3.2](https://github.com/gtonkinhill/panaroo)
- [Sourmash v4.8.4](https://github.com/sourmash-bio/sourmash)

```bash
# Generate checkM results. Change the expected extension of the assemblies accordingly in the -x flag
checkm lineage_wf -t ${cpus} -x fa --tab_table genomes_dir checkm_output > checkm_out.tsv
echo "genome,completeness,contamination,strain_heterogeneity" > checkm_quality.csv
cat checkm_out.tsv | cut -f1,11,12,13,14 | sed 's/\t0\t/\t/;/Completeness/d;s/\t/,/g;s/,/.fa,/' >> checkm_quality.csv

# Create the representative species clusters
dRep dereplicate -g genomes_dir/* \
   -p ${cpus} \
   -pa 0.9 -sa 0.95 -nc 0.30 \
   -cm larger -comp 50 -con 5 \
   --genomeInfo checkm_quality.csv \
   drep_output

# Taxonomic labelling of the representative genomes
gtdbtk classify_wf \
   --cpus ${cpus} \
   --pplacer_cpus ${pp.cpus} \
   --genome_dir drep_output/dereplicated_genomes/ \
   --extension fa \
   --skip_ani_screen \
   --out_dir gtdbtk_results
cat gtdbtk_summary_bac120 gtdbtk_summary_arc53 > gtdbtk.summary.tsv

# Functional annotation of all genomes according to the dRep table drep_output/data_tables/Cdb.csv
# Consider that some genomes in genomes_dir could be discarded due to checkM QC
cat Cdb.csv | cut -d',' -f1 | sed '/genome/d;s/\..*//' > catalogue_genomes.list 
for genome in $(cat catalogue_genomes.list); do (
   prokka genomes_dir/$genome.fa \
      --cpus ${cpus} \
      --kingdom 'Bacteria' \
      --outdir $genome\_prokka \
      --prefix $genome \
      --force \
      --locustag $genome
); done

# Second step of functional annotation on amino acid sequences
emapper.py -i $genome\_prokka/$genome.faa \
   --database eggnog/data/eggnog.db \
   --dmnd_db eggnog/data/eggnog_proteins.dmnd \
   --data_dir eggnog/data/ \
   -m diamond \
   --no_file_comments \
   --cpu ${cpus} \
   --dbmem \
   -o $genome\_eggnog

# Computing the pangenome per species cluster
# Preparing inputs to run Panaroo. Move the gff files from prokka_results to a new directory called `gff_files`
mkdir gff_files
for genome in $(cat catalogue_genomes.list); do (
   mv prokka_results/$genome\_prokka/$genome.gff gff_files
); done

# Creating the directory structure to run Panaroo
panaroo_inputs_builder_custom.py \
   --drep_clstrs drep_output/data_tables/Cdb.csv \
   --gff_files ./gff_files \
   --derep_genomes drep_output/dereplicated_genomes/

# Running Panaroo on each species cluster
panaroo -i ${DIR}/*.gff \
    -t ${cpus} \
    -o ${DIR}_panaroo \
    --clean-mode strict \
    --merge_paralogs \
    --core_threshold 0.90 \
    --threshold 0.90 \
    --family_threshold 0.5 \
    --no_clean_edges

# Generating the core_genes tables
for dir in $(ls -d rep_*_panaroo | sed 's/rep_//;s/_panaroo//'); do ( echo "Processing "$dir && cd rep_$dir\_panaroo/ && get_core_genes.py -i gene_presence_absence.Rtab -o core_genes.txt ); done

# Tidying up the directories
mkdir panaroo_results
mv *_panaroo panaroo_results
rm -r rep_*

# Computing the pangenome tables for the BioSIFTR pipeline
pangenomeDB_builder_custom.py \
   --drep_clstrs drep_output/data_tables/Cdb.csv \
   --derep_genomes drep_output/dereplicated_genomes/ \
   --pfam_dat Pfam-A.hmm.dat.gz \
   --panaroo_path panaroo_results \
   --prokka_path prokka_results

mkdir functional_profiles_DB
mv *_clstr.tsv functional_profiles_DB

# Compute the KEGG modules' completeness tables as explained above

# In addition, you will need to create the indexed database for Sourmash kmer = 21 of your representative genomes and (optionally) the bwamem2 indexed database
sourmash sketch dna -p scaled=1000,k=21 --name "$prefix" $path -o $prefix.fna.sig
sourmash index sourmash_species_representatives_k21.sbt.zip *sig

cat drep_output/dereplicated_genomes/* > bwa_reps.fna
bwa-mem2 index bwa_reps.fna

```


### 4. Functional annotation benchmark
The optimisation parameters were set up in the MGnify BioSIFTR pipeline, and the tool was used to generate the functional profiles to benchmark against the ground truth. Results were generated using `-core_mode true` and `--core_mode false`.

```bash
# Running the BioSIFTR pipeline on the synthetic communities
nextflow run ebi-metagenomics/biosiftr --input samplesheet.csv --outdir pan_shallow_results --biome chicken-gut-v1-0-1 --run_bwa true 

nextflow run ebi-metagenomics/biosiftr --input samplesheet.csv --outdir core_shallow_results --biome chicken-gut-v1-0-1 --run_bwa true --core_mode true

# Integrating functional annotation matrices. The following commands were run to integrate Pfam, KOs and KEGG modules annotation.
matrix_integrator.py --input ground_truth/kegg_counts.tsv pan_shallow_results/integrated_annotation/bwa_kos_matrix.tsv core_shallow_results/integrated_annotation/bwa_kos_matrix.tsv pan_shallow_results/integrated_annotation/sm_kos_matrix.tsv core_shallow_results/integrated_annotation/sm_kos_matrix.tsv --output integrated_kos.txt

# Compute performance metrics from functions matrices
matrix2performance.py --input integrated_kos.txt --output metrics_kos.txt

```

The concatenated table of performance generated in this section is available in the [data/optimisation
/function/](https://github.com/EBI-Metagenomics/shallow_shotgun_paper/tree/main/data/optimisation/function) directory of this repo.



<a name="sec2"></a>
## Section 2. Pipeline validation on real data
The BioSIFTR pipeline was challenged on three different gut biomes of real data: human ([PRJDB11444](https://www.ebi.ac.uk/ena/browser/view/PRJDB11444)), mouse ([PRJEB74255](https://www.ebi.ac.uk/ena/browser/view/PRJEB74255)), and junglefowl ([PRJEB46806](https://www.ebi.ac.uk/ena/browser/view/PRJEB46806)). 
Results were compared with other strategies for functional prediction:

- 16S rRNA amplicon
   - [Picrust2](https://github.com/picrust/picrust2)
   - [MicFunPred](https://github.com/microDM/MicFunPred)
- Deep shotgun
   - Assembly with [SPAdes](https://github.com/ablab/spades) or [Megahit](https://github.com/voutcn/megahit) and functional annotation with [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper)
- Shallow-shotgun (artificially subsampled from deep-shotgun)
   - [SHOGUN](https://github.com/knights-lab/SHOGUN)
   - [BioSIFTR pipeline](https://github.com/EBI-Metagenomics/biosiftr)

<p align="center" width="100%">
   <img src="visuals/validation.png" width="100%"/>
</p>


### 1. Functional annotation and comparative metagenomics analysis
The raw reads of 16S rRNA amplicon were processed using QIIME to generate ASVs and PICRUSt2 and MicFunPred were used to generate the functional profiles.
This is an example for one of the three data sets (Red junglefowl; accession PRJEB46806). The code used for the two data sets is identical, except for the accession numbers (PRJDB11444, PRJEB74255) and the options "--p-trunc-len-f" and "--p-trunc-len-r" in the denoising step.
Briefly, PRJEB46806 reads were trimmed to 230 bp on both fw and rv reads (as in this example), PRJDB11444 reads were trimmed to 280 bp on the fw read and 220 bp on the rv read, and PRJEB74255 reads were trimmed to 245 bp on the fw read and 220 bp on the rv read.

```bash
# Code used to generate ASVs, taxonomic labelling, and function inference from PRJEB46806 amplicon data
#The filereport "filereport_read_run_PRJEB46806_tsv.txt" should be downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJEB46806?show=reads
grep "MiSeq" "filereport_read_run_PRJEB46806_tsv.txt" | awk '{print $11}' | awk -v RS=';' 1 > PRJEB46806_files.txt
wget -i PRJEB46806_files.txt
ls *.fq.gz | awk -F '_' '{print $1}' | sort | uniq | awk -v pwd="$PWD" -v c="/" -v OFS="\t" '{print $1, pwd c $1 "_16S_R1.fq.gz", pwd c $1 "_16S_R2.fq.gz"}' > PRJEB46806_manifest.txt
sed -i '1s/^/sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n/' PRJEB46806_manifest.txt

#all files are phred-33 coded

#Import data as a Qiime2 artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path PRJEB46806_manifest.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

#Check sequence quality for trimming
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

#Denoise to ASVs with dada2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 230 \
  --p-trunc-len-r 230 \
  --p-min-overlap 20 \
  --p-pooling-method pseudo \
  --p-n-threads 8 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

#Export ASV-table to tab-delimited format
qiime tools export \
  --input-path table.qza \
  --output-path qiime_exports

biom convert -i qiime_exports/feature-table.biom -o qiime_exports/ASV-table.tsv --to-tsv

#Export representative sequences
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path qiime_exports

#Predict functionality with MicFunPred and Picrust2
MicFunPred_run_pipeline.py -i qiime_exports/ASV-table.tsv -r qiime_exports/dna-sequences.fasta -o micfunpred_output -t 8 --contrib -v
picrust2_pipeline.py -i qiime_exports/feature-table.biom -s qiime_exports/dna-sequences.fasta -o picrust2_output -p 8

#Export files
tar -cvf PRJEB46806_KO_output.tar.gz picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz micfunpred_output/KO_metagenome/KO_metagenome.tsv.gz

# Code used to generate ASVs, taxonomic labelling, and function inference from PRJEB46806 amplicon data

```

Deep-shotgun raw reads were filtered by quality, decontaminated, and assembled before functional annotation.

```bash
# Deep-shotgun raw-reads quality control (fastp)
fastp --in1 sample_1.fq.gz --in2 sample_2.fq.gz --out1 sample_1.filt.fq.gz --out2 sample_2.filt.fq.gz --json sample.fastp.json --html sample.fastp.html --thread 6 --detect_adapter_for_pe

# Deep-shotgun decontamination of human (GCF_000001405.40), phyX (NC_001422.1) and chicken (GCF_016699485.2)
bwa-mem2 mem -M -t 16 hg38.fa sample_1.filt.fq.gz sample_2.filt.fq.gz | \
samtools view -@ 16 -f 4 -F 256 -uS - | \
samtools sort -@ 16 -O bam - -o sample_sorted.bam
samtools bam2fq -@ 8 -1 sample_1.hg38.fq.gz -2 sample_2.hg38.fq.gz -0 /dev/null -s /dev/null -n sample_sorted.bam

bwa-mem2 mem -M -t 16 phiX174.fna sample_1.hg38.fq.gz sample_2.hg38.fq.gz | \
samtools view -@ 16 -f 4 -F 256 -uS - | \
samtools sort -@ 16 -O bam - -o sample_sorted.bam
samtools bam2fq -@ 8 -1 sample_1.phix.fq.gz -2 sample_2.phix.fq.gz -0 /dev/null -s /dev/null -n sample_sorted.bam

bwa-mem2 mem -M -t 16 chicken.fna sample_1.phix.fq.gz sample_2.phix.fq.gz | \
samtools view -@ 16 -f 4 -F 256 -uS - | \
samtools sort -@ 16 -O bam - -o sample_sorted.bam
samtools bam2fq -@ 8 -1 sample_1.decont.fq.gz -2 sample_2.decont.fq.gz -0 /dev/null -s /dev/null -n sample_sorted.bam

# HQ and decontaminated reads assembly with SPAdes or MEGAHIT when memory exceeds 200 G
spades.py --only-assembler --meta --threads 16 --memory 120 -1 sample_1.decont.fq.gz -2 sample_2.decont.fq.gz --continue -o sample
megahit -1 sample_1.decont.fq.gz -2 sample_1.decont.fq.gz -t 16 --out-dir sample

# Filtering out contigs length < 500 bp
assembly_filter.py contigs.fasta 500

# EggNOG annotation on contigs of length > 500 bp
emapper.py -i 500_contigs.fasta --itype metagenome --database eggnog/data/eggnog.db --dmnd_db eggnog/data/eggnog_proteins.dmnd --data_dir eggnog/data/ -m diamond --no_file_comments --cpu 16 --dbmem

# Extracting KOs annotation (emapper evalue threshold = 1E-10)
eggnog2KOs.py --eggnog sample.emapper.annotations --sample sample

```

We generated shallow-shotgun datasets from the deep-shotgun HQ decontaminated samples having >10 M raw reads by random subsampling using the [seqtk](https://github.com/lh3/seqtk) tool. Then we ran the BioSIFTR pipeline and the SHOGUN tool to generate taxonomic and functional profiles. A [SHOGUN docker image](https://quay.io/repository/microbiome-informatics/shogun_knights_lab_1.0.8) developed in this work is available on quay.io.

```bash
# Generating artificial shallow-shotgun data (sequencing yields = 100000, 500000, 1000000, 1500000, 2000000)
SEED=$(shuf -i 0-999 -n 1)
seqtk sample -s $SEED sample_1.decont.fq.gz 100000 > sample_100000_1.fq
seqtk sample -s $SEED sample_2.decont.fq.gz 100000 > sample_100000_2.fq

# Running the BioSIFTR pipeline
nextflow run ebi-metagenomics/biosiftr --input samplesheet.csv --outdir junglefowl_results --biome chicken-gut-v1-0-1 --run_bwa true

# Running SHOGUN. Requires transforming fastq to fasta, fixing read names and concatenating the reads
fq2fa.py sample_100000_1.fq.gz sample_100000_R1
cat sample_100000_R1.fasta | sed 's/_/sub/;s/_/|/;s/_//g;s/|/_/' > sample_100000.fasta 
cat sample_100000_R2.fasta | sed 's/_/sub/;s/_/|/;s/_//g;s/|/_/' >> sample_100000.fasta 

shogun pipeline --input sample_100000.fasta --database shogun/databases/ --output sample_100000 --level species --function -t 16

```

Functional tables generated by the different methods were integrated, and a comparative metagenomics analysis was performed in R using the [ordination_plots_shallow.ipynb](https://github.com/EBI-Metagenomics/shallow_shotgun_paper/tree/main/notebooks) notebook available in this repository.

```bash
matrix_integrator.py --input picrust2.tsv micfunpred.tsv deep.tsv shallow.tsv shogun.tsv --output mixed_vals_kos.tsv

# Transforming into absence/presence matrix
relab2abspres.py --matrix mixed_vals_kos.tsv --output presabs_kos.tsv
```

KEGG Orthologues annotation tables generated by each of the different methods are available in the [data/validation](https://github.com/EBI-Metagenomics/shallow_shotgun_paper/tree/main/data/validation) directory of this repo.


### 2. Taxonomic annotation

This section corresponds to the post-processing of taxonomic tables generated by the different methods. When necessary, the transformation from NCBI taxonomy to GTDB was made to GTDB r207.

Taxonomic profiles were generated from the deep-shotgun HQ decontaminated reads using the [mOTUs pipeline](https://github.com/EBI-Metagenomics/motus_pipeline) and taxonomy was transformed into GTDB names using the strategy suggested by [Alessio Milanese](https://github.com/motu-tool/mOTUs/wiki/GTDB-taxonomy-for-the-mOTUs).

```bash
# Amplicon data
#Classify rep sequences against Greengenes2
wget http://ftp.microbio.me/greengenes_release/2024.09/2024.09.taxonomy.asv.nwk.qza
wget http://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.fna.qza

qiime greengenes2 non-v4-16s \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --i-backbone 2024.09.backbone.full-length.fna.qza \
  --o-mapped-table table_gg2.qza \
  --o-representatives rep-seqs_gg2.qza

qiime greengenes2 taxonomy-from-table \
  --i-reference-taxonomy 2024.09.taxonomy.asv.nwk.qza \
  --i-table table_gg2.qza \
  --o-classification table_gg2_taxonomy.qza

#Summarize and export files
qiime tools export \
  --input-path table_gg2.qza \
  --output-path qiime_exports_gg2

biom convert -i qiime_exports_gg2/feature-table.biom -o qiime_exports_gg2/ASV-table_gg2.tsv --to-tsv

qiime tools export \
  --input-path rep-seqs_gg2.qza \
  --output-path qiime_exports_gg2

qiime tools export \
  --input-path table_gg2_taxonomy.qza \
  --output-path qiime_exports_gg2

#Export files
tar -czvf PRJEB46806_gg2_output.tar.gz qiime_exports_gg2/

# Processing ASVs table to add clean taxonomy labels. Removing GG2 (_[0-9]+) monophyletic identifiers in taxonomic ranks and aggregating names
asv2taxo.py --asv_table qiime_exports_gg2/ASV-table_gg2.tsv --names_table taxonomy.tsv 

# Transform amplicon count tables into relative abundance, removing singletons. Discarded singletons and doubletons: 16
counts2relab.py --count_table amplicon_taxo.tsv --output relab_amplicon_taxo.tsv

# Deep-shotgun data. Taxonomic annotation with the mOTUs pipeline
nextflow run motus_pipeline/main.nf --mode paired --paired_end_forward sample_1.decont.fq.gz --paired_end_reverse sample_2.decont.fq.gz --sample_name sample 

# Transform mOTUs NCBI labels into GTDB taxonomy
motus2gtdb.py --input sample_merged.fastq.tsv --mapping mOTUs_3.0.0_GTDB_tax.tsv --sample sample

# Transform the mOTUs count tables into relative abundance, removing singletons. Discarded singletons and doubletons: 127
counts2relab.py --count_table motus_gtdb.tsv --output relab_motus_taxo.tsv

# Shallow-shotgun data. From BioSIFTR results
# Aggregating shallow_mapping results at species level
shallow2aggr.py --taxonomy_table sm_taxo_matrix.tsv
shallow2aggr.py --taxonomy_table bwa_taxo_matrix.tsv
matrix_integrator.py --input sm_taxo_matrix.tsv bwa_taxo_matrix.tsv --output shallow_taxo.tsv

# Shallow-shotgun data. From SHOGUN results
# Transforming labels from NCBI to GTDB r207. GTDB metadata files were downloaded from https://data.gtdb.ecogenomic.org/releases/release207/207.0/. Discarded singletons and doubletons: 807
shogun2gtdb.py --shogun_tax shogun_taxo.tsv --metadata_r207 bac120_metadata_r207.tsv.gz ar53_metadata_r207.tsv.gz --output shogun_gtdb.tsv
counts2relab.py --count_table shogun_gtdb.tsv --output relab_shogun_gtdb.tsv


# Integrating all the taxonomic tables. All use GTDB r207 (Greengenes2 2022.10, motus)
matrix_integrator.py --input relab_amplicon_taxo.tsv relab_motus_taxo.tsv shallow_taxo.tsv relab_shogun_gtdb.tsv --output taxonomy_relab.tsv

# Aggregating to different taxonomic ranks
ranks_aggregator.py --input taxonomy_relab.tsv --output taxonomy_relab

# Transforming into presence/absence
relab2abspres.py --matrix taxonomy_relab_domain.tsv --output taxonomy_presabs_domain.tsv

```

Ordination plots and other comparative metagenomics analysis were computed in R using the code in the [ordination_plots_biosiftr.ipynb](https://github.com/EBI-Metagenomics/biosiftr_extended_methods/tree/main/notebooks) notebook available in this repository.

