# Intro

This is a work-in-progress repository

## Workflow

This is a general workflow for ancestral reconstruction with links/tags to more detailed sections

1. **[Genome_selection](./Notes/Genome_selection.md)**  First we select the genomes to use and download from NCBI/collect. This includes checking completeness [checkm2_results](./Results/checkm2%20results.xlsx) and average nucleotide identity  [fastani_results](./Results/fastani_results_filtered.txt)
 #completeness #ani

2. **[Annotation](Notes/Annotation.md)** Genomes here annotated using [MetaErg2](https://github.com/kinestetika/MetaErg)  
#metaerg #annotation

3. [Species Tree](Notes/species_tree.md)

    1. Tree_of_mags
    2. Run raxml (a bit time-consuming)
    3. Reroot the species tree

4. Gene trees

    1. Finding Orthologs
    2. Run MSA
    3. IQtree fast bootstrap to produce [gene trees](Notes/Gene_Trees.md)
    4. Convert everything into ALEobjects

5. [Gene-tree-species-tree Reconciliation](./Notes/Gene-tree-species-tree-reconcile.md)