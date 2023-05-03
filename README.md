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
        a. tree_of_mags uses conserved single copy marker gene HMMs - same as used gtdb ~ 50. e.g, rRNA, etc
    2. Run raxml (a bit time-consuming)
    3. Reroot the species tree

4. [Gene trees](./Notes/Gene_Trees.md)

    1. Finding [Orthologs](./Notes/Gene_Trees.md#run-orthologs) 
        a. this is essentially a homology search. This script is going to produce 1000s of clusters of orthologs and paralogs
    2. Run MSA [using clustalo](./Notes/Gene_Trees.md#create-msas)
    3. IQtree fast bootstrap to produce [gene trees](Notes/Gene_Trees.md)
    4. Convert everything into ALEobjects

5. [Gene-tree-species-tree Reconciliation](./Notes/Gene-tree-species-tree-reconcile.md)