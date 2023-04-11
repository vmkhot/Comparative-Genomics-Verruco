# Analyses after the reconciliation

After figuring out which root is the most likely in [gene-tree-species-tree reconciliation](./Gene-tree-species-tree-reconcile.md#which-root-is-correct), we can run more analyses to interpret the data. 

From section 3.9: We can evaluate the robustness of the root

## Robustness Checks

From [the textbook](https://link.springer.com/protocol/10.1007/978-1-0716-2691-7_9)

>In traditional phylogenetics, site stripping—the removal of fast-evolving [49] or compositionally biased [50, 51] sites—is often
used to assess the robustness of phylogenetic signal. The rationale
is that the information in these sites can be misleading (or at least
difficult to model adequately), for example, favoring long branch or
compositional attraction over real evolutionary signal. Relationships in the overall tree that derive their support from fast-evolving
or compositionally biased sites are then considered less reliable than
relationships that are robust to the removal of these sites. An
analogous situation exists in the context of root or topology inference at the phylogenomic level: gene families differ in terms of size,
evolutionary rate, and rates of DTL events, and some families may
be easier to model than others given the methods available...  
By analogy to the traditional case, relationships that derive their support from particular
subsets of the data (e.g., gene families with very high loss or transfer
rates) might be considered less reliable than those which obtain
support across the full distribution of families.

We can test the robustness of each root by evaluating the impact of high rates of duplication, transfer or loss on each root using [DTL_ratio_analysis_ML_diff.py](https://github.com/ak-andromeda/ALE_methods/blob/main/DTL_ratio_analysis_ML_diff.py) script from [ALE_methods](https://github.com/ak-andromeda/ALE_methods)

**You need 2 other files for this script to run [species_list_demo.txt](../species_list_demo.txt) (listing all the species) and a [roots_to_test.txt](../roots_to_test.txt) (paths to the roots directories)**

```bash
# create a conda environment for the packages
conda create -n ALE

# install necessary packages
conda install -c conda-forge pyparsing pandas=1.5.3 seaborn
```

Run the script

```bash
#!/bin/bash
#SBATCH --job-name=other          # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=10           # Number of CPU cores per task
#SBATCH --mem=10G                   # Job memory request
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=other%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ale

# test each condition separately
# root 3 given here is the "most likely root"

python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 LS
python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 L
python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 DS
python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 D
python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 TS
python3 ../ALE_methods/DTL_ratio_analysis_ML_diff.py root3 T
```

This script produces 3 plots. Not sure which to look at yet?? 2 plots are the difference in the summed likelihood from the "chosen root".  
[T_110423_difference](../Results/Root_robustness_checks/T_110423_difference.png) is all the gene familes  
[T_110423_high_species_difference](../Results/Root_robustness_checks/T_110423_difference.png) is only gene families with > 50% species representation  

[T_110423_summed](../Results/Root_robustness_checks/T_110423_summed.png) What is this one mean?

E.g. Transfer rates  
Gene families are *sequentially* removed based on their tranfer rates. High transfer rates are removed first. And the summed difference in likelihood (ML_diff)
from the maximum likelihood root (dotted line at 0) recalculated at each iteration.  

**x axis**: e.g. "T_rank" is the ML_diff at each iteration of removal of one gene family. That's why it goes from 0 to 7294 (total # of gene families in our analyses)

The rest of ML difference plots are gathered together [here (all gene families)](../Results/Root_robustness_checks/All_difference.png) and [here (species rep > 50%)](../Results/Root_robustness_checks/All_high_species_difference.png)

*What does this finally mean for my chosen root? Is it not robust enough to the removal of high rates of DTL/DS,TS,LS?*