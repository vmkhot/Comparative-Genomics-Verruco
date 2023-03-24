# Gene Trees

Using all 47 selected genomes [genome_selection](./Genome_selection.md) and their amino acid files produced by metaerg2.

## Run Orthologues script

### Download CompGenomics package

Refer to [installation instructions](./species_tree.md#download-compgenomics-package)

### Run Orthologs

The default minimum representation for clusters is 3 genes. This can be changed with the `--minimum_representation 0` argument  
#MCL #Clustering #trees #gene_tree

```bash
# Make your output directory
mkdir ./orthologs_out
```

```bash
#!/bin/bash
#SBATCH --job-name=orthologs           
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                   
#SBATCH --mem=50G                   
#SBATCH --time=24:00:00              
#SBATCH --output=orthologs%j.log      
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate CompGenomics

orthologues --input_dir ../metaerg/all_genomes_faa/ --output_dir ./orthologs_out --cpus 30 --minimum_representation 3
```

This resulted in 9000+ clusters with >=3 representatives in each cluster.

Visually inspect random clusters to see if the genes grouped together have similar annotations.

## Create MSAs

Yes, we will need to run over 9000+ multiple sequence alignments (woo!)  

```bash
#!/bin/bash
#SBATCH --job-name=clustalo           
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=30           
#SBATCH --mem=50G                   
#SBATCH --time=48:00:00
#SBATCH --output=clustalo%j.log       

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

#clustal omega is in the CompGenomics conda environment
conda activate CompGenomics

# run MSA on each fasta file as a loop
for fn in ./orthologs_out/*.faa;
do
    base="${fn:16:-4}"
    echo $base
    clustalo -i $fn -o clustalo_out/${base}_MSA.fa -v --threads 30
done
```

IQtree will error out if the alignments have less than 3 sequences in them so have to check this:

```bash
# count headers in each fasta using grep
grep -c ">" ../orthologs/clustalo_out/ > num_genes.list

# sort to see them
sort -n -t ':' -k2 num_genes.list | less

# move them
mv ../orthologs/clustalo_out/101_MSA.fa  ../orthologs/clustalo_out/small_clusters/
```

## Run IQtree

#iqtree #substitutionmodel #trees  
The following command does approximately this: take each gene cluster > find the best substitution model > start with a random tree with best fit model > generate this 10,000 times (stored in .ufboot) > the type of tree that is generated the most often is a consensus tree (stored in .contree) > take the log likelihood of this consensus tree > produce maximum likelihood tree (stored in .treefile)

```bash
#!/bin/bash
#SBATCH --job-name=iqtree           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30           # Number of CPU cores per task
#SBATCH --mem=100G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=cpu2023,cpu2019         # Partition to bigmem because lots of memory required
#SBATCH --output=iqtree%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate iqtree

#modelfinder
#iqtree2 -s ../orthologs/clustalo_out/1300_MSA.fa -m MF -redo


# need to run iqtree as a loop. giving it the whole directory doesn't work in "-s".
# I created and ram 9 slurm scripts. There is also an MPI option but requires a diff compilation

# -MFP: run model finder and use that to make tree
# -madd LG+C20,LG+C60: mixed models, LG is general amino acid subs model
# -B 10000: 10,000 bootstraps
# -T, --mem: thread and memory constraints on the program. have to give memory option, otherwise it will time-out
for fn in ../orthologs/clustalo_out/1*_MSA.fa;
do
    iqtree2 -s $fn -m MFP -madd LG+C20,LG+C60 -B 10000 -wbtl -T 40 --mem 100G
done
```

After running this program, you will get the following files in your alignments folder.  
|File|Description|
|-|-|
0_MSA.fa|original alignment
0_MSA.fa.bionj|rapid neighbour joining tree
0_MSA.fa.ckp.gz|
0_MSA.fa.contree| consensus tree
0_MSA.fa.iqtree| iqtree report
0_MSA.fa.log| log file
0_MSA.fa.mldist| maximum likelihood distances
0_MSA.fa.model.gz| model fitting results
0_MSA.fa.splits.nex|
0_MSA.fa.treefile| maximum likelihood tree (best consensus maximum likelihood gene tree)
0_MSA.fa.ufboot| 10,000 bootstrapped trees (to be used later)

Next we use the bootstrapped trees (.ufboot) and the [species tree](./species_tree.md) to do the [gene tree-species tree reconciliation](./Gene-tree-species-tree-reconcile.md)

