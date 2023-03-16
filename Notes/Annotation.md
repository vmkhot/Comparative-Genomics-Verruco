# Annotation

Genomes were annotated using MetaErg2  
#annotation #metaerg

[To download MetaErg using apptainer](https://github.com/vmkhot/Metagenome-workflows/blob/main/Illumina-Short-Reads/Annotation.md#metaerg-20)

To run on slurm

```bash
#!/bin/bash
#SBATCH --job-name=metaerg2_1      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=40            # Number of CPU cores per task
#SBATCH --mem=100G                    # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=metaerg2_1run%j.log     # Standard output and error log
pwd; hostname; date



apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --database_dir /work/ebg_lab/referenceDatabases/metaerg_database/ --path_to_signalp /work/ebg_lab/referenceDatabases/metaerg_database/ --path_to_tmhmm /work/ebg_lab/referenceDatabases/metaerg_database/ --contig_file ../genomes2run/genomes_1 --rename_contigs --cpus 40 --file_extension .fna --force antismash --force aragorn --output_dir ./genomes_1_output
```
