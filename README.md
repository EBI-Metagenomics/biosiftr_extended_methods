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
  2. Taxonomic prediction power
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

[InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) was used to generate each of the synthetic communities from [MetaChic](https://entrepot.recherche.data.gouv.fr/dataset.xhtml?persistentId=doi:10.15454/FHPJH5) chicken-gut genomes catalogue using representative genomes only. The genome contigs were renamed to make easier to track back the synthetic genes origin after mapping

```bash
# Renaming the contigs of all representatives
assembly_rename.py rep_genome.fa rep_genome



```





