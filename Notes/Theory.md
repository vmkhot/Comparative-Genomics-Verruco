# Theory

All the theory we learnt in class

## Week 1

Reading?
What is the problem?

## Week 2

You're comparing proteins from different genomes

what is a protein in a bioinformatics sense?

What is a small protein?
50aa (small) -- 300aa (avg) -- 1000aa(large)

A very small bacterium: ~1000um average size
the size of membrane width: ~10nm - this translates to 50aa  

Protein sequences in biofinformatics are long words of amino acids

Proteins in reality have 3D structure, an active site and a shape (particular fold)

Some amino acids are more important than others in the sequence to maintain the same activity, folding and structure (i.e maintaining function).

So proteins can be upto 80% different and still share the same function

**Homologs:** (BROAD) two proteins that share some function. e.g. an iron transporter and calcium transporter

**Orthologs:** (IN SAME SPECIES) Two proteins in same a common ancestor

**Paralogs:** An homologous protein that has arisen in organism B but is not shared with the ancestor. A second copy of the gene that comes out because of HGT or duplication

### Practice

|A|B|
|---|---|
|2000|2500|

How many homologues exist?
This is a trick question! because there might be more than 1. So we are actually looking for families of homologs

|homolog H1|
|---|
|A53|
|A100|
|A10004|
|B93|

So the fact that we find 3 homologues of H1 in organism A is suspect in our application. Our application is to figure out if the organism has gained or lost in genes in its adaptation.

**To find homologues** we use BLAST using a concatenated database of organism A proteins and organism B proteins.

```bash
cat A.fa B.fa > catAB.fa
makeblastdb
blastp catAB.fa catAB.fa
```

Can also use DIAMOND for scaling up

To cluster proteins, use MCL? (markov clustering). MCL clusters by the graph that is produced by blastp

Can also use usearch or CDhit but we have to define the identity threshold for clustering - but what to use? - this is problematic

MCL has a network expansion parameter? has a range from 1-5 - have to run multiple times to figure out how tight you want the clustering to be?

|clustering results|||
|---|---|---
|1|(very narrow clustering) 4 clusters of 3 proteins|
|2|4 clusters of 3 proteins
|3|4 clusters of 3 proteins
|4|1 clusters of 7 proteins
|5|(very broad clustering)1 clusters of 7 proteins

How many organisms are represented?
You want the most number of organisms represented and the fewest paralogs.  
*How do we know which are paralogs?*  
We take the protein with the highest identity to the cluster? this is a crude approach.  
Really we need to make a gene tree and species tree for ancestral reconstruction to figure out which proteins are orthologs and which were gained and lost

**Questions for Marc:** Frame shift correction brought the completion to 95%. GTDB closest relative completion is 89%.??
Did FastANI between my genome and all other from Opitutae (1400+) genomes and no results. Closest genome ANI from gtdbtk was to GCA_003566375.1 at 77.01%

*Can I use genomes which are not complete?* depends how incomplete they are. 5% is likely okay but 30% incomplete will be missing many genes. #completeness

## Week 3


We need to create a heirarchy of cluters as the inner clusters have different levels of identity to each other.
e.g. 50%: A-Mg, B-Mg, A-Ca, B-Ca --> 1 cluster --> this is the parent cluster --> **1 ortholog and 2 paralogs**
    70%: (A-Mg, B-Mg), (A-Ca, B-Ca) --> 2 clusters --> children of the 50% --> **2 orthologs**
    90&: (A-Mg, B-Mg), A-Ca, B-Ca --> 1 cluster and 2 singletons --> children of the 70% --> **1 ortholog**

so the most ideal cluster for the above transporters is at 70%

Need to carefully choose which genomes to cluster with target genome

[Genome_selection](./Genome_selection.md)

## Week 4

What's next?
what's missing?

how to make a species tree?
tree_of_mags script:

- recruit single-copy marker genes with HMMs (genes) - HMMs are same as GTDBtk
- redo the orthologs with these single copy marker genes (homologs)
- align each cluster (multiple MSA)
- concatenate all MSAs (1 MSA)
- inspect MSA visually
- We do this: Use the MSAs to create the tree - RAXML, fasttree, 
- Inspect the tree

Run CheckM2 to see the completeness of all the genomes/MAGs we are going to use.

 #MCL #Clustering [Trial](./Trial_w_2genomes.md)

MCL can also cluster gene-expression data (transcriptomes and proteomes)

MCL output interpretation
it starts with the largest cluster and then goes to the smaller clusters. 

Ideally, we want clusters with only 1 gene from each organism so for two organisms - clusters of 2 genes.

How to break apart clusters of 200+? These could be families of orthologs that are closely related or actual paralogs. --> We run mcl multiple times with different inflation values and then pick the best inflation value for each gene family. There's no good inflation value for all the genes because soem genes evolve faster than others.

how to pick the best clusters?
priority for picking:

1. most orthologs
2. least paralogs

what are the next steps? After picking the best clusters, you will have groups of orthologs. You can pick out the groups representing single-copy-marker-genes and use this gene content to make a robust evolutionary species tree. 

Then you can ask similar questions about other genes by constructing gene trees from the ortholog groups

### Chapter 9 from Environmental Microbial Evolution
 #trees #rootingtrees

***Rooting Species Trees Using Gene Tree-Species Tree
Reconciliation*** by Brogan J. Harris, Paul O. Sheridan, Adria´n A. Davı´n, Ce´cile Gubry-Rangin, Gergely J. Szo¨llo˝si, and Tom A. Williams

Each branch has a branch length.
Each tree has a scale bar - this represents the # of substitution per site. (e.g. 0.01 - means 1 substitution per site for 100 sequences) - this is what determines the evolutionary distance. (e.g. length of 0.2 - 20 substitutions per 100 sequences - this is very fast evolution)

#### ALE

 #ALE  
Amalgamated Likelihood Estimation  
Dating trees - how does this happen?

Molecular information + fossils to establish when something happens. e.g. origin of photosynthesis
Read Chapter 9 - they use IQtree

Molecular clock hypothesis: the rate of molecular evolution (sequence substitution) is constant over time. This is obviously not strictly true as 1 clade might evolve faster than another. How: bacterial doubling time and mutation?

Strategies for rooting:

- Selection of outgroup: a branch from the last common ancester of subject clade. If outgroup is very distantly related - problematic tree artifacts occur, such as long branch attraction
- Midpoint root: the midpoint of the longest tip-to-tip distance. Follows strict molecular clock
- Relaxed Molecular Clock (RMC): models evolutionary rate variation and incorporates fossil calibrations, but is computationally intensive
- non-stationary, non-reversible: not commonly used

Reconciliation methods, on the other hand, aim to extract
rooting information from phylogenetic discord between gene trees and the species tree, incorporating signal from gene duplications and transfers. Gene ree-species tree reconciliation works by drawing gene trees into the species tree using a series of gene and species-level events, such as gene duplication,transfer, and losses. Parsimony-based reconciliation methods seek reconciliations that imply the fewest events given prespecified event costs, while probabilistic methods rely on a generative model of the evolution of gene trees along the species tree.

## Week 5

- Need to download metaerg 2.0
- Run genomes/MAGs through Metaerg2
- update tree-of-mags script
- Use .faa files to run tree-of-mags
- inspect MSA
- download raxml
- create tree Raxml
- inspect tree --> compare to expectation from other papers or gtdb
- reroot the tree on the outgroup. B/c trees are always unrooted - the program makes no assumptions on who evolved first.
- Need to rename all the genomes in the tree (need script for this)

### Annotation

 #annotation #metaerg  
 Where does annotation go wrong?

#### Gene prediction

People are interested in protein-coding genes and therefore use Prodigal for gene prediction. However, prodigal does not think about other genomic features such as CRISPRs, tRNAs, rRNA genes, repeat sequences

Therefore, we have to predict all the other features first and **MASK** the features for prodigal by replacing the region with Ns (NNNN)

#### Functional and Taxonomic annotation

We start with protein coding genes then using homology we predict their function. Annotators use both

- Homology: BLAST or HMMs
- Signal peptides - 20aa that exists at the N-terminus that tells the ribosome that this needs to go outside the cell or where to go
- Transmembrane helix - membrane bound proteins

The problem is with HMMs and BLAST. Proteins that are the same are called slightly different things which makes it difficult to interpret as a human reader. Every dataabase has a problem with naming based on the knowledge bank of the annotators or who built the database.

We also want to know the taxonomy PER GENE b/c it can point us to where the gene might have originated from.
Metaerg: GTDB taxonomy AND NCBI annotations + custom HMMs

Metaerg is run on the DNA files

## Week 6

- run CheckM2 on all genomes
- ALE protocol - amalgamated likelihood expectation? protocol. From Environmental Microbial Evolution

### Substitution models

#substitutionmodel  
For IQtree, have to choose a substitution model, which affects the final tree significantly.  
**closely related sequences:**
only the codons are changing  
**distantly related sequences:**
lots of things are changing so have to look specifically at how similar the amino acids are.

What does Laura Hug do for species tree? PROTGAMMALG for RAXML

For next week:

- Look at model choices:
  - try a couple to see difference
  - Hug et al. uses PROTGAMMALG
  - look at ancestral reconstruction papers
- Final species tree
- Run orthologues script using metaerg output
  - check visually 10 of the ortholog groups to see if the annotations match for all the genes in a group
- Run alignments using different aligners
  - clustalO, mafft, muscle
- IQtree fastbootstrap - how long does it take? use the time option on the command line?

Things to do:

1. Run metaerg again on all folders
2. update tree_of_mags and orthologues again
3. Run species tree again with PROTGAMMALG in RAXML
4. Inspect tree
5. Find out about genomes in the sister clade (DAMON??)
6. Run Orthologues script with --minimum_representation argument

## Week 8

Check if the command "parallel" can be used with IQTREE?

What are the next steps? Convert all the ufboots to ALEobjects.  
Rooting:  
Next we try to root the species tree using likelihood (step 3.7 in the textbook). How? Set the root at some hypothetical locations and the ALE script estimates the likelihood of each root.  
How?: by drawing genes trees into the species tree - which one is correct? one method is parsimony - the least number of duplication and transfer events. Another is probabilistic methods - what is this? Reading for next week.

inferring gene content:  
inferring gene content at each position in the tree {gene loss, gain, duplication and transfer}