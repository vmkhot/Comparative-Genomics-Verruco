# Gene tree species tree reconciliation

First have to download alesuite. I downloaded the latest docker image and use apptainer to access

```bash
apptainer pull docker://boussau/alesuite:latest
```
AleObserve to convert the bootstrapped trees (.ufboot)  to ale objects 

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