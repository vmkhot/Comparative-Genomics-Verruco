# Intro

This is a work-in-progress repository

## Workflow

This is a general workflow for ancestral reconstruction with links/tags to more detailed sections

1. **[Genome_selection](./Notes/Genome_selection.md)**  First we select the genomes to use and download from NCBI/collect. This includes checking completeness [checkm2_results](./Results/checkm2%20results.xlsx) and average nucleotide identity  [fastani_results](./Results/fastani_results_filtered.txt)
 #completeness #ani

2. **[Annotation](./Notes/Annotation.md)** Genomes here annotated using [MetaErg2](https://github.com/kinestetika/MetaErg)  
#metaerg #annotation

3. [Species Tree](Notes/species_tree.md)

4. Finding Orthologs

5. IQtree fast bootstrap to produce [gene trees](Notes/Gene_Trees.md)

6. Convert everything into ALEobjects