# Species Tree
#trees #species_tree #rootingtrees
Using all 47 selected genomes [genome_selection](./Genome_selection.md) and their amino acid files produced by metaerg2.

## Run tree_of_mags script

### Download CompGenomics package

```bash
conda create -n CompGenomics

conda activate

# close the git repository
git clone https://github.com/kinestetika/ComparativeGenomics.git
```

go into `ComparativeGenomics/src/comparative_genomics/tree_of_mags.py` file and change `java -jar /bio/bin/BMGE/src/BMGE.jar` to just `bmge`. bmge was installed using conda into the CompGenomics environment

```bash
# continue installation
cd ComparativeGenomics/

python -m build

# change version number (0.7) below accordingly
pip install --upgrade dist/comparative_genomics-0.13.tar.gz
```

To update, have to delete the `ComparativeGenomics` directory and redo installation from the git clone.

### Run tree_of_mags

```bash
#!/bin/bash
#SBATCH --job-name=tree_of_mags           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30           # Number of CPU cores per task
#SBATCH --mem=100G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem          # Partition to bigmem because lots of memory required
#SBATCH --output=tree_of_mags%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate CompGenomics

# tree of mags command
# --keep_identical
# --min_representation: The minimum number of taxa represented in a gene for it to be kept (default =0)
# --min_frequency: The minimum fraction of genes a taxon should be represented in for the taxon to be kept. (default =0)

tree_of_mags --dir ../metaerg/all_genomes_faa/ --file_extension ".faa" --cpus 30 --keep_identical
```

This produces a file in the "alignments" directory called `alignments/concatenated_alignment` - then this concatenated alignment was used to run raxml-ng to produce the species tree. the model `--model LG+G` parameter should be same as `PROTGAMMALG`. `LG` is the default for amino acid alignments and `G` is the `GAMMA`  
#RAXML #species_tree

### Make the species tree

```bash
#!/bin/bash
#SBATCH --job-name=raxml
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --output=raxml%j.log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate raxml
#parse argument will just check the alignment and estimate the computation time/resources required
# --threads auto{35} :optimize parallelization within 35 threads (what was requested from slurm)

raxml-ng --msa alignments/concatenated_alignment  --threads auto{35} --model LG+G


#old command
#raxml-ng --msa alignments/concatenated_alignment.fasta  --threads 8 --model LG+G+F --tree rand{2}

pwd; hostname; date
```

Resulting tree is [species_tree](../Results/tree2_concatenated_alignment.raxml.bestTree)

### Reroot the species tree
#rootingtrees

For the ALE reconciliation, we require >=3 rooted trees with different potential roots. To produce these there are 2 options:  

1. Manually choose root on [itol](https://itol.embl.de/) and save as newick
2. Use a python script [root_tree_in_position.py](https://github.com/ak-andromeda/ALE_methods/blob/main/root_tree_in_position.py)

For this, I manually chose the roots in ITOL based on metadata about the genomes:  

1. Outgroup chosen earlier [newick_tree_outgroup_1](../Results/rerooting_trees/newick_tree_outgroup_1.txt)
2. Midpoint root (chosen by itol) [newick_tree_midpoint_2](../Results/rerooting_trees/newick_tree_midpoint_2.txt)
3. Clade of Soda lake genomes first (included the Opitutales genome) [newick_tree_sodalake_3](../Results/rerooting_trees/newick_tree_sodalake_3.txt)
4. Clade of Freshwater genomes first [newick_tree_freshwater_4](../Results/rerooting_trees/newick_tree_freshwater_4.txt)
5. Clade of Marine genomes first [newick_tree_marine_5](../Results/rerooting_trees/newick_tree_marine_5.txt)