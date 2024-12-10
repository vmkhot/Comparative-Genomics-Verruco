Number of proteins:

```bash
[vmkhot@arc blastp]$ grep -c ">" *.faa
GCA_003566375.1_proteins.faa:3212
Opitutales_own_genome_proteins.faa:5249
concat_own_GCA_003566375.1_proteins.faa:8461
```

### [Clustering with MCL](http://micans.org/mcl/#:~:text=The%20MCL%20algorithm%20is%20short,in%20bioinformatics%20and%20other%20disciplines.)

 #MCL #Clustering

MCL creates an undirected network.
How does the algorithm work? [Read here](https://towardsdatascience.com/markov-clustering-algorithm-577168dad475)

change the blast format to pairwise connections with query, subject and evalue

```bash
cut -f 1,2,11 seq.cblast > seq.abc
```

```bash
#!/bin/bash
#SBATCH --job-name=mcl           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20           # Number of CPU cores per task
#SBATCH --mem=50G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --out=mclclustering%j.log       # Standard index and error log

pwd; hostname; date

# mcxload to create network files
# the evalues are converted to a value between 1 and 200 with a negative log10 transformation
mcxload -abc ownGenome_GCA_003566375.1.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o ownGenome_GCA_003566375.1.mci -write-tab ownGenome_GCA_003566375.1.tab

#run mcl clustering with different inflation parameters
mcl ownGenome_GCA_003566375.1.mci -I 1
mcl ownGenome_GCA_003566375.1.mci -I 1.2
mcl ownGenome_GCA_003566375.1.mci -I 1.4
mcl ownGenome_GCA_003566375.1.mci -I 2
mcl ownGenome_GCA_003566375.1.mci -I 4
mcl ownGenome_GCA_003566375.1.mci -I 6


##label the index_out.uts using mcxdump
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I10 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I10
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I12 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I12
 mcxdump -icl out.ownGenome3_GCA_003566375.mci.I14 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I14
 mcxdump -icl out.ownGenome3_GCA_003566375.mci.I20 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I20
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I30 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I30
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I40 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I40
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I50 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I50
mcxdump -icl out.ownGenome3_GCA_003566375.mci.I60 -tabr ownGenome3_GCA_003566375.tab -o dump.ownGenome3_GCA_003566375.mci.I60

#alternatively, you can generate tab files directly from mcl
#mcl seq.mci -I 1.4  -use-tab ownGenome_GCA_003566375.1.tab
#mcl seq.mci -I 2  -use-tab ownGenome_GCA_003566375.1.tab
#mcl seq.mci -I 4  -use-tab ownGenome_GCA_003566375.1.tab
#mcl seq.mci -I 6  -use-tab ownGenome_GCA_003566375.1.tab
```

We could use more than 4 inflation values to get more granularity. Here I used 8 just to see how the clusters would change.

```bash
I=1 [mcl] 45 clusters found
I=1.2 [mcl] 1304 clusters found
I=1.4 [mcl] 1674 clusters found
I=2 [mcl] 1990 clusters found
I=3 [mcl] 2227 clusters found
I=4 [mcl] 2368 clusters found
I=5 [mcl] 2443 clusters found
I=6 [mcl] 2482 clusters found
```

_How low does the % identity have to be to reach similar number of clusters as mcl?_  
Only clusters with two or more should be counted

### [Clustering with vsearch instead](https://github.com/torognes/vsearch)

 #Clustering  
 Again, we cluster at multiple sequence identities to find the best fit of clusters. E.g. for --id 0.8, vsearch will cluster all sequences that are at least 80% identical.

```bash
# download with Conda
conda install -c bioconda vsearch

source ~/miniconda3/etc/profile.d/conda.sh
conda activate vsearch

# clusters.uc : file describing which sequences are clustered together. 
vsearch --cluster_fast genome3_GCA_003566375_genes.fna --centroids vsearch_proteins_centroids_0.8.faa --id 0.8 -uc clusters_0.8.uc

vsearch --cluster_fast genome3_GCA_003566375_genes.fna --centroids vsearch_proteins_centroids_0.9.faa --id 0.9 -uc clusters_0.9.uc

vsearch --cluster_fast genome3_GCA_003566375_genes.fna --centroids vsearch_proteins_centroids_0.95.faa --id 0.95 -uc clusters_0.95.uc

vsearch --cluster_fast genome3_GCA_003566375_genes.fna --centroids vsearch_proteins_centroids_0.98.faa --id 0.98 -uc clusters_0.98.uc
```

Results

```bash
# number of clusters at each identity
16976 clusters_0.8.uc
17050 clusters_0.9.uc
17085 clusters_0.95.uc
17137 clusters_0.98.uc
```

## Comparative Genomics Scripts

[Scripts from Marc](https://github.com/kinestetika/ComparativeGenomics)

For the Comparative Genomics scripts you require a python environment with python > 3.10

```bash
conda create -n CompGenomics python=3.10.0

conda activate CompGenomics
```

Installation

```bash
# go to directory you want to install in
cd ~/Programs/
# download from github
git clone https://github.com/kinestetika/ComparativeGenomics.git

#build the package
cd ComparativeGenomics
python -m build

# install the package (check version and adjust)
pip install --upgrade dist/comparative_genomics-0.4.tar.gz
```

Required programs for running the scripts:
Prodigal
hmmer suite
clustal omega
BMGE

To run the orthologues script

```bash
#!/bin/bash
#SBATCH --job-name=orthologs           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30           # Number of CPU cores per task
#SBATCH --mem=100G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem                      # Partition to bigmem because lots of memory required
#SBATCH --output=orthologs%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate CompGenomics

orthologues --input_dir ./all_genomes --output_dir ./orthologs_out --cpus 30
```
