# Gene tree species tree reconciliation

First have to download alesuite. I downloaded the latest docker image and use apptainer to access

```bash
apptainer pull docker://boussau/alesuite:latest
```

[AleObserve](https://github.com/ssolo/ALE) to convert the bootstrapped trees (.ufboot)  to ale objects 

```bash
#!/bin/bash
#SBATCH --job-name=aleO          # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=40           # Number of CPU cores per task
#SBATCH --mem=50G                   # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=aleO_1_%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

# used GNU parallel to run multiple instances at once. 

ls ../../orthologs/clustalo_out/UFBOOT/1*.ufboot| parallel apptainer exec -B /work/ebg_lab/eb/ ~/data/Programs/alesuite_latest.sif ALEobserve {} burnin=1000
```

## Reconciliation

#ALE #rootingtrees

Now we combine the gene trees (.ale) objects and the [rooted species trees](./species_tree.md#reroot-the-species-tree) using the `ALEml_undated` command from alesuite.

The reconciliation for each root needs to be run in its own directory. The following script needs to be run in each root directory

separators='|' : describes genome|gene  
fraction_missing= : describes missing part of each genome based on checkm2 [in sheet2](../Results/checkm2%20results.xlsx)

```bash
#!/bin/bash
#SBATCH --job-name=ml_undated          # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=40           # Number of CPU cores per task
#SBATCH --mem=50G                   # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=ml_undated_%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

# try on small gene tree first
apptainer exec -B /work/ebg_lab/eb/ ~/data/Programs/alesuite_latest.sif ALEml_undated newick_tree_outgroup_1.txt ../../orthologs/clustalo_out/UFBOOT/100_MSA.fa.ufboot.ale separators='|' fraction_missing='../fraction_missing.txt'


# using parallel to run all the ale objects
parallel "apptainer exec -B /work/ebg_lab/eb/ ~/data/Programs/alesuite_latest.sif ALEml_undated newick_tree_outgroup_1.txt {} separators='|' fraction_missing='../fraction_missing.txt'" ::: ../../orthologs/clustalo_out/UFBOOT/new_folder/*.ale
```

Three files are produced for each gene. Example results are [here](../Results/ml_undated_example/)

## Which root is correct?

To identify which root is the most likely, we construct a likelihoods table using the following script from [ALEmethods](https://github.com/ak-andromeda/ALE_methods). 

```bash
python3 ../ALE_methods/write_consel_file_p3.py root1 root2 root3 root4 root5 > likelihoods_table.mt
```
Then we use this program CONSEL to consolidate the likelihoods table into human-readable-format. 

To download [CONSEL](http://stat.sys.i.kyoto-u.ac.jp/prog/consel/)

```bash
# get directory clone from github
git clone https://github.com/shimo-lab/consel.git

cd consel/src
make
mkdir ../bin
cp catass catci catmt catpv catrep consel makerep makermt randrep treeass ../bin/
```

To run CONSEL

```bash
makermt likelihoods_table.mt
consel likelihoods_table
catpv likelihoods_table > au_test_out
```

Results are in [au_test_out](../Results/au_test_out.txt) - what do these mean?  
The roots are ranked by likelihood and p-value(au column). In this case  
p < 0.05: throw it out  
p > 0.05: cannot be rejected  
*This means that the null hypothesis is that the root cannot be rejected*

In our case, the highest ranked au is **root3: sodalake first**  
 So perhaps the organism transitioned from soda lake to marine and not the other way around. Midpoint root and Marine root could also not be rejected but freshwater and outgroup roots (both p < 0.05) are out.
