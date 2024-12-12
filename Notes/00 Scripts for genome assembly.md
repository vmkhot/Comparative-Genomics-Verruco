## QC
### Porechop
```bash
#!/bin/bash
#SBATCH --job-name=porechop     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=porechop_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate porechop

for file in ../no_sample/20220708_1153_MN26797_FAQ98508_b39622cb/fastq_pass/*.fastq.gz;
do
    sample="${file:0:-9}"
    echo $sample
    porechop -i $file -o ${sample}_trimmed.fastq.gz -t 16
done
```

### nanofilt
```bash
#!/bin/bash
#SBATCH --job-name=nanofilt     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --output=nanofilt_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate nanofilt 

#gunzip -c all_reads_trimmed.fastq.gz | NanoFilt -q 10 | gzip > all_reads_trimmed_HQ.fastq.gz

NanoPlot -t 30 --fastq all_reads_trimmed.fastq.gz --plots dot --legacy hex
```

## Assembly

### Canu
```bash
#!/bin/bash
#SBATCH --job-name=canu     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --output=canu_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

#canu is a pipeline that also does read correction so it requires untrimmed, raw reads
#conda activate canu

canu -useGrid=remote -p opitutales -d opitutales genomeSize=2.3m maxInputCoverage=100 -nanopore ../all_passed_reads.fastq.gz
```
### Minimap2+miniasm
```bash
#!/bin/bash
#SBATCH --job-name=miniasm     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --output=miniasm_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate minisuite

minimap2 -x ava-ont ../all_reads_trimmed_HQ.fastq.gz ../all_reads_trimmed_HQ.fastq.gz | gzip > ./minimap.paf.gz

miniasm -f ../all_reads_trimmed_HQ.fastq.gz ./minimap.paf.gz > miniasm.gfa

awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fasta

assembly-stats ./miniasm.fasta
```
### Raven

```bash
#!/bin/bash
#SBATCH --job-name=raven     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --output=stats_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate raven

#raven -t 30 -p 4 ../all_reads_trimmed_HQ.fastq.gz

conda activate minisuite
assembly-stats ./raven.fasta
```

### Trycycler
```bash
#!/bin/bash
#SBATCH --job-name=trycycler     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node   
#SBATCH --ntasks=1                   # Run a single task        
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --time=30:00:00              # Time limit hrs:min:sec
#SBATCH --output=trycycler_%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate trycycler

#to cluster the assemblies

trycycler cluster --assemblies ../*.fasta --reads ../../all_reads_trimmed_HQ.fastq.gz --out_dir ./output

#to reconcile the clusters
trycycler reconcile --reads ../../all_reads_trimmed_HQ.fastq.gz --cluster_dir output/cluster_001

trycycler partition --reads ../../all_reads_trimmed_HQ.fastq.gz --cluster_dirs output/cluster_001

trycycler consensus --cluster_dir output/cluster_001
```

## Polishing
```bash
#!/bin/bash
#SBATCH --job-name=polish     # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=polish%j.log     # Standard output and error log
pwd; hostname; date


#MAPPING
source ~/miniconda3/etc/profile.d/conda.sh

#conda activate minisuite

#polish1
#long-read mapping
#minimap2 -ax map-ont ../02_Assembly/7_final_consensus.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_consensus_assembly_minimap2.sam

#conda activate racon
#racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_consensus_assembly_minimap2.sam ../02_Assembly/7_final_consensus.fasta  > racon1.fasta

#polish2
conda activate minisuite
minimap2 -ax map-ont ./racon1.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon1_minimap2.sam

conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon1_minimap2.sam racon1.fasta  > racon2.fasta

#polish3
conda activate minisuite
minimap2 -ax map-ont ./racon2.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon2_minimap2.sam

conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon2_minimap2.sam racon2.fasta  > racon3.fasta

#polish4
conda activate minisuite
minimap2 -ax map-ont ./racon3.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon3_minimap2.sam

conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon3_minimap2.sam racon3.fasta  > racon4.fasta
```

## CheckM2

```bash
#!/bin/bash
#SBATCH --job-name=checkm      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=100G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=checkm%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

#conda activate checkm2

checkm2 predict -t 30 -x fasta --input ./ --output-directory ./CheckM2
```

## Taxonomy
```bash
#!/bin/bash
#SBATCH --job-name=phylo      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=phylo_test%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

#conda activate gtdbtk

gtdbtk classify_wf -x fasta --cpus 20 --genome_dir ../Consensus_Genome/ --out_dir ./ --pplacer_cpus 20 --write_single_copy_genes --keep_intermediates
```

