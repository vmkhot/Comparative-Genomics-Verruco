# [[01 Genome Selection]]

## selecting genomes from the gtdb metadata file

```python
#! /home/vmkhot/miniconda3/envs/new-working-env/bin/python3.12


import pandas as pd
from pandas import DataFrame
import numpy as np
import pyarrow

gtdb_df = pd.read_csv("./bac120_metadata_r220.tsv.gz", sep="\t", compression='gzip',engine='pyarrow')
#print(gtdb_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER GTDB BY CRITERIA AND SORT SO THAT BEST QUALITY GENOMES ARE ON TOP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

verruco_df = gtdb_df.loc[gtdb_df['gtdb_taxonomy'].str.contains('p__Verrucomicrobiota;') &
                                      (gtdb_df['checkm_completeness'] >= 90) & (gtdb_df['checkm_contamination'] <= 5) & 
                                      (gtdb_df['checkm2_completeness'] >= 90) & (gtdb_df['checkm2_contamination'] <= 5) &
                                      (gtdb_df['checkm_strain_heterogeneity'] == 0)]



# split gtdb taxonomy
verruco_df[['domain','phylum','class','order','family','genus','species']]=verruco_df['gtdb_taxonomy'].str.split(';', expand = True)

# sort the columns to select genus reps
verruco_df.sort_values(by=['genus','gtdb_type_species_of_genus','mimag_high_quality','checkm2_completeness','checkm2_contamination'],
                                  ascending=[True,False,False,False,True], inplace=True)
print(verruco_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUBSETTING DATAFRAMES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select 1 class level rep for each class (except c__Verrucomicrobiae)
verruco_class_reps_df = verruco_df.loc[~verruco_df['class'].str.contains('c__Verrucomicrobiae')].groupby('class').first().reset_index()
print(verruco_class_reps_df)

# select 2 order level rep for each order (except o__Opitutales)
verruco_order_reps_df = verruco_df.loc[verruco_df['class'].str.contains('c__Verrucomicrobiae') & ~verruco_df['order'].str.contains('o__Opitutales')].groupby('order').head(2).reset_index()
print(verruco_order_reps_df)

# select 3 family level rep for each family (except f__T3Sed10-336)
# need them to be unique genomes (to be informative for phylogeny)
verruco_family_reps_df_best = verruco_df.drop_duplicates('gtdb_taxonomy', keep='first') \
.loc[verruco_df['class'].str.contains('c__Verrucomicrobiae') & verruco_df['order'].str.contains('o__Opitutales') & ~verruco_df['family'].str.contains('f__T3Sed10-336')] \
.groupby('family').head(3).reset_index()

print(verruco_family_reps_df_best)

verruco_family_reps_df_genus = verruco_df.drop_duplicates('genus', keep='first') \
.loc[verruco_df['class'].str.contains('c__Verrucomicrobiae') & verruco_df['order'].str.contains('o__Opitutales') & ~verruco_df['family'].str.contains('f__T3Sed10-336')] \
.groupby(['family']).head(3).reset_index()
# .groupby(['family'], as_index=False, group_keys=False,
#          ).apply(
#              lambda x: x.sample(min(3,len(x)))
#          ).reset_index()

print(verruco_family_reps_df_genus)

# concatenate in order of unique genera and then best representative. That way in the group, there's unique genera first in each group. If not 3 unique then takes the next best representative for the family
df_merge_verrcu_family = pd.concat([verruco_family_reps_df_genus,verruco_family_reps_df_best]).drop_duplicates('species') \
.groupby('family',sort=False).head(3) \
.sort_values(by=['family','genus']) \
.drop(['index'],axis=1).reset_index()

print(df_merge_verrcu_family)

# select all members of T3Sed10-336 family
# some are duplicates
verruco_T3S_df = verruco_df.loc[verruco_df['family'].str.contains('f__T3Sed10-336')]
print(verruco_T3S_df)


# select 5 random refseq genomes from Plancto phylum
plancto_df = gtdb_df.loc[gtdb_df['gtdb_taxonomy'].str.contains('p__Planctomycetota;') & 
                          gtdb_df['accession'].str.contains('GCF') &
                          (gtdb_df['checkm_completeness'] > 98) & (gtdb_df['checkm_contamination'] <= 2) & 
                          (gtdb_df['checkm_strain_heterogeneity'] == 0)]

plancto_df[['domain','phylum','class','order','family','genus','species']]=plancto_df['gtdb_taxonomy'].str.split(';', expand = True)
plancto_df.sort_values(by=['genus','gtdb_type_species_of_genus','mimag_high_quality','checkm_completeness','checkm_contamination'],
                                  ascending=[True,False,False,False,True], inplace=True)

gtdb_plancto_genus_df = plancto_df.groupby('genus').first().reset_index()
# print(gtdb_plancto_genus_df)

# select 5 random genus level reps with a reproducible seed of 123
plancto_sample_df = gtdb_plancto_genus_df.sample(n=5, random_state=123)
print(plancto_sample_df)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONCATENATE DATAFRAMES AND CLEAN
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pd.concat will concatenate unordered columns also by matching them as keys
all_selected_genomes_gtdb = pd.concat([verruco_class_reps_df,verruco_order_reps_df,df_merge_verrcu_family,verruco_T3S_df,plancto_sample_df],ignore_index=True, sort=True).drop_duplicates(subset='accession')

# fix up accessions
all_selected_genomes_gtdb[['RS_GB','accession']]=all_selected_genomes_gtdb['accession'].str.split('_', n=1, expand = True)

#replace all GCF with GCA
all_selected_genomes_gtdb['accession'] = all_selected_genomes_gtdb['accession'].str.replace('GCF','GCA')   

print(all_selected_genomes_gtdb)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WRITE OUT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_selected_genomes_gtdb.to_csv('final_genomes_selection_Verruco.tsv',header=True,sep='\t',index=False)

download_genomes_df = all_selected_genomes_gtdb[['accession','ncbi_assembly_name']].copy()
download_genomes_df.to_csv('genomes_to_download.tsv',index=False, header=False)



```

Add the soda lake genomes in manually

## downloading genomes from ncbi (genbank)

```python
#! /home/vmkhot/miniconda3/envs/new-working-env/bin/python3.12


#! usr/bin/python3

# genbank ftp link format
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/644/195/GCA_001644195.2_ASM164419v2/GCA_001644195.2_ASM164419v2_genomic.fna.gz

"""
1. Get genbank accessions and assembly names
2. Split into variables
3. Create ftp link
4. Download to a genomes dir
"""

from pathlib import Path
import os, shutil               # for directory and file control commands like cp and mv
import urllib.request

Path("./fna_files").mkdir(exist_ok=True) # have to make the output directory if it doesn't exist. same behaviour as (mkdir -p)

ftp_path = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"

with open ("genomes_to_download.csv", 'r') as accessions:
    for line in accessions:
        cols=line.strip('\n').split(',')
        genbank = cols[0]
        assembly = cols[1].str.replace(' ','_')
        first = genbank[4:7]
        second = genbank[7:10]
        third = genbank[10:13]
        #print(first, second, third)
        ftp_retrieve_genomic = str.join('/',(ftp_path, str(first), str(second), str(third),f"{genbank}_{assembly}",f"{genbank}_{assembly}_genomic.fna.gz"))
        ftp_retrieve_genomic = ftp_retrieve_genomic.replace(' ','_')
        print(ftp_retrieve_genomic)
        try:
            urllib.request.urlretrieve(ftp_retrieve_genomic,f"./fna_files/{genbank}_{assembly}_genomic.fna.gz")
        except:
            print(genbank + "...not found fna")
            pass


# check
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/376/875/GCA_002376875.1_ASM237687v1/GCA_002376875.1_ASM237687v1_genomic.fna.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/376/875/GCA_002376875.1_ASM237687v1/GCA_002376875.1_ASM237687v1_genomic.fna.gz


```

# [02 Annotation, Orthologs, MSA and Gene Trees]()
## metaerg

metaerg in comparative analysis mode should give genome annotations + identify orthogroups 
```bash

#!/bin/bash
#SBATCH --job-name=metaerg2      
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=40            
#SBATCH --mem=100G                    
#SBATCH --time=48:00:00              
#SBATCH --output=metaerg2_%j.log     
#SBATCH --mail-user=vmkhot@ucalgary.ca     
#SBATCH --mail-type=END                         

pwd; hostname; date


apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --database_dir /work/ebg_lab/referenceDatabases/metaerg_db_V214/ --path_to_signalp /work/ebg_lab/referenceDatabases/metaerg_db_V214/ --path_to_tmhmm /work/ebg_lab/referenceDatabases/metaerg_db_V214/ --contig_file /work/ebg_lab/eb/Varada/Verruco_redo/download_genomes/fna_files --rename_genomes --cpus 40 --file_extension .fna --mode comparative_genomics --output_dir /work/ebg_lab/eb/Varada/Verruco_redo/metaerg_annotated_genomes/annotations/

```

### split CDS files into directories with number of sequences

```python
#! /home/vmkhot/miniconda3/envs/new-working-env/bin/python3.12

import glob, os
from pathlib import Path
import shutil

#for i in range(1,44):
for file in Path(f"./untrimmed_MSA").glob("*.faa"):
    num = len([1 for line in open(file.resolve()) if line.startswith(">")])
    #print(num)
    if num < 4 :
        new_path = str(file).replace(f"untrimmed_MSA", "seq_lt_4_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 4 <= num < 25 :
        new_path = str(file).replace(f"untrimmed_MSA", "seq_lt_25_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 25 <= num < 50 :
        new_path = str(file).replace(f"untrimmed_MSA", "seq_25_to_50_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 50 <= num < 75:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_50_to_75_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 75 <= num < 100:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_75_to_100_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 100 <= num < 200:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_100_to_200_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 200 <= num < 300:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_200_to_300_trim")    # make the directory first
        shutil.move(file, new_path)
    elif 300 <= num < 500:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_300_to_500_trim")    # make the directory first
        shutil.move(file, new_path)
    else:
        new_path = str(file).replace(f"untrimmed_MSA", "seq_mt_500_trim")    # make the directory first
        shutil.move(file, new_path)

```

Only work with more than 4 sequences from now forward

Split files in the directories into subdirs

```bash
i=0; for f in seq_lt_25/*.faa; do d=seq_lt_25/dir_$(printf $((i/70+1))); mkdir -p $d; mv "$f" $d; let i++; done

i=0; for f in seq_25_to_50/*.faa; do d=seq_25_to_50/dir_$(printf $((i/30+1))); mkdir -p $d; mv "$f" $d; let i++; done

i=0; for f in seq_50_to_75/*.faa; do d=seq_50_to_75/dir_$(printf $((i/30+1))); mkdir -p $d; mv "$f" $d; let i++; done

i=0; for f in seq_mt_75/*.faa; do d=seq_mt_75/dir_$(printf $((i/20+1))); mkdir -p $d; mv "$f" $d; let i++; done

```


## MSA with famsa

Ran famsa in a loop 

famsa_job_template.sh
```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=<TIME>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<MEM>
#SBATCH --output=<OUTPUT_FILE>


echo "FAMSA 2"

source /vast/ri65fin/miniconda3/etc/profile.d/conda.sh

conda activate famsa  # Activate conda environment 

mkdir <OUT_DIR>

for f in <IN_DIR>/*.faa
do
    newname=$(basename $f .faa)
    famsa -t 0 $f <OUT_DIR>/${newname}_MSA.faa
done

```

famsa_submit_jobs.sh
```bash
#!/bin/bash

# Loop over the job indices
for i in {1..102}
do
    # Replace placeholders in the job template with actual values
    job_name="famsa_$i"
    input_dir="../30_CDS/seq_lt_25/dir_$i"
    output_file="famsa_out_$i.log"
    output_dir="./seq_lt_25/dir_$i"
    job_file="famsa_seq_lt_25_$i.sbatch"
    
    sed -e "s|<JOB_NAME>|$job_name|g" \
        -e "s|<OUTPUT_FILE>|$output_file|g" \
        -e "s|<IN_DIR>|$input_dir|g" \
        -e "s|<OUT_DIR>|$output_dir|g" \
        -e "s|<TIME>|48:00:00|g" \
        -e "s|<NODES>|1|g" \
        -e "s|<NTASKS>|1|g" \
        -e "s|<CPUS>|10|g" \
        -e "s|<MEM>|10G|g" \
        30_famsa.sbatch > $job_file

    # Submit the job
    sbatch $job_file
done
```

Famsa job template

```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=<TIME>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<MEM>
#SBATCH --output=<OUTPUT_FILE>
#SBATCH --mail-user=vmkhot@ucalgary.ca
#SBATCH --mail-type=END


echo "FAMSA 2"

source /home/vmkhot/miniconda3/etc/profile.d/conda.sh

conda activate famsa  # Activate conda environment 

mkdir -p <OUT_DIR>

for f in <IN_DIR>/*.faa
do
    newname=$(basename $f .faa)
    famsa -t 0 $f <OUT_DIR>/${newname}_MSA.faa
done

```


## trim alignments using trimAL

trimal_job_template.sh
```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=<TIME>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<MEM>
#SBATCH --output=<OUTPUT_FILE>


echo "trimAl"

mkdir <OUT_DIR>

for f in <IN_DIR>/*.faa
do
    newname=$(basename $f .faa)
    trimal -in $f -out ${newname}_trimmed.faa -gappyout
done

```

trimal_submit_jobs.sh (submit multiple times for the folders)
```bash
#!/bin/bash

# Loop over the job indices
for i in {1..1}
do
    # Replace placeholders in the job template with actual values
    job_name="trimal_$i"
    input_dir="./seq_mt_75/dir_$i"
    output_file="trimal_seq_mt_75_$i.log"
    output_dir="./seq_mt_75_trimmed/dir_$i"
    job_file="trimal_seq_mt_75_$i.sbatch"
    
    sed -e "s|<JOB_NAME>|$job_name|g" \
        -e "s|<OUTPUT_FILE>|$output_file|g" \
        -e "s|<IN_DIR>|$input_dir|g" \
        -e "s|<OUT_DIR>|$output_dir|g" \
        -e "s|<TIME>|24:00:00|g" \
        -e "s|<NODES>|1|g" \
        -e "s|<NTASKS>|1|g" \
        -e "s|<CPUS>|20|g" \
        -e "s|<MEM>|10G|g" \
        20_trimal.sbatch > $job_file

    # Submit the job
    sbatch $job_file
done
```

trimal job template

```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=<TIME>
#SBATCH --nodes=<NODES>
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<MEM>
#SBATCH --output=<OUTPUT_FILE>
#SBATCH --mail-user=vmkhot@ucalgary.ca
#SBATCH --mail-type=END

echo "trimAl"

mkdir -p <OUT_DIR>

for f in <IN_DIR>/*.faa
do
    newname=$(basename $f .faa)
    trimal -in $f -out <OUT_DIR>/${newname}_trimmed.faa -gappyout
done

```

## Gene trees using iqtree2

Now run iqtree2 on the trimmed MSAs (this will take ages)
We will run everything with LG+C20+G model

iqtree_submit_jobs_25_50.sh (submit multiple, one for each directory)
```bash
#!/bin/bash

mkdir -p seq_25_to_50

# Loop over the job indices
# Loop over the job indices
for i in {1..102}
do
    # Replace placeholders in the job template with actual values
    job_name="iqtree_$i"
    input_dir="../40_MSA/seq_lt_25_trimmed/dir_$i"
    output_dir="./seq_lt_25/dir_$i"
    output_file="./sbatch_scripts_logs/iqtree_seq_lt_25_%j.log"
    job_file="./sbatch_scripts_logs/iqtree_job_seq_lt_25_$i.sbatch"
    
    sed -e "s|<JOB_NAME>|$job_name|g" \
        -e "s|<OUTPUT_FILE>|$output_file|g" \
        -e "s|<IN_DIR>|$input_dir|g" \
        -e "s|<OUT_DIR>|$output_dir|g" \
        -e "s|<TIME>|120:00:00|g" \
        -e "s|<NTASKS>|1|g" \
        -e "s|<CPUS>|5|g" \
        -e "s|<MEM>|2G|g" \
        -e "s|<GUIDE>|LG+F+G|g" \
        -e "s|<MODEL>|LG+C20+G|g" \
        10_iqtree.sbatch > $job_file

    # Submit the job
    sbatch $job_file
done
```

iqtree job template

```bash
#!/bin/bash
#SBATCH --job-name=<JOB_NAME>
#SBATCH --time=<TIME>
#SBATCH --ntasks=<NTASKS>
#SBATCH --cpus-per-task=<CPUS>
#SBATCH --mem=<MEM>
#SBATCH --output=<OUTPUT_FILE>

echo "iqtree 2"

source /home/vmkhot/miniconda3/etc/profile.d/conda.sh

conda activate iqtree  # Activate conda environment 
 
mkdir -p <OUT_DIR>
mkdir -p <OUT_DIR>/guide_trees

for f in <IN_DIR>/*_trimmed.faa
do
    newname=$(basename $f .faa)
    iqtree2 -s $f -m <GUIDE> -T <CPUS> --mem <MEM> -pre <OUT_DIR>/guide_trees/pmsf_guide_${newname}

    iqtree2 -s $f -m <MODEL> -ft <OUT_DIR>/guide_trees/pmsf_guide_${newname}.treefile -mwopt  -B 5000 -wbtl -T <CPUS> --mem <MEM> -pre <OUT_DIR>/pmsf_${newname}
done


# some trees were remade using the option -keep-ident in the second command to keep less than 4 non-unique sequences in the analysis
```

# 03 Species tree

```bash
#!/bin/bash
#SBATCH --job-name=species_tree      
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=40            
#SBATCH --mem=100G                    
#SBATCH --time=24:00:00              
#SBATCH --output=species_tree_%j.log     
#SBATCH --mail-user=vmkhot@ucalgary.ca     
#SBATCH --mail-type=END

# make gtdb-tk alignment using single-copy markers

source /home/vmkhot/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk-2.4.0

gtdbtk identify --genome_dir ../20_metaerg_annotated_genomes/annotations/fna/ --out_dir ./11_gtdbtk_results_renamed/identify --extension fna --cpus 40

gtdbtk align --identify_dir ./11_gtdbtk_results_renamed/identify  --out_dir ./11_gtdbtk_results_renamed/align  --cpus 40

conda deactivate

# IQTREE

conda activate iqtree  # Activate conda environment 

mkdir -p ./20_iqtree/

# unzip the MSA
gunzip ./11_gtdbtk_results_renamed/align/align/gtdbtk.bac120.user_msa.fasta.gz

# run iqtree2
iqtree2 -T 40 -m TEST -s ./11_gtdbtk_results_renamed/align/align/gtdbtk.bac120.user_msa.fasta -B 1000 -alrt 1000 --keep-ident -pre ./20_iqtree/species_tree_renamed

```

# 04 ALE Reconciliation

### AleObserve to create ALE objects

Ran a burn-in of 500 only - because only have 5000 gene trees in the distributions

```bash
#!/bin/bash
#SBATCH --job-name=aleobserve      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=50G                    
#SBATCH --time=48:00:00              
#SBATCH --output=aleobserve_%j.log
#SBATCH --error=aleobserve_%j.err   
#SBATCH --mail-user=vmkhot@ucalgary.ca     
#SBATCH --mail-type=END


pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

# used GNU parallel to run multiple instances at once. 

find ../50_gene_trees/seq_* -name *.ufboot | parallel apptainer exec -B /work/ebg_lab/eb/ ~/data/Programs/alesuite_latest.sif ALEobserve {} burnin=500
```

###  Ale_ml_undated

```bash
#!/bin/bash
#SBATCH --job-name=ml_undated      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150G                    
#SBATCH --time=48:00:00              
#SBATCH --output=ml_undated_%j.log
#SBATCH --error=ml_undated_%j.err   
#SBATCH --mail-user=vmkhot@ucalgary.ca     
#SBATCH --mail-type=END


pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

# used GNU parallel to run multiple instances at once. 

parallel "apptainer exec -B /work/ebg_lab/eb/ ~/data/Programs/alesuite_latest.sif ALEml_undated species_tree_renamed_rooted.tree {} fraction_missing=fraction_missing.txt" ::: /work/ebg_lab/eb/Varada/Verruco_redo/50_gene_trees/ale_objects/*.ale
```

# 05 Further analyses

get_gene_events_all_nodes.py

```python
## This script will get gene family events at all nodes (including the leaves) and output both with and without the cutoff

import os
import sys
import glob
import pandas as pd
from pathlib import Path

def get_copies_at_node(rec_file):

    """ Function to get the copies of the gene family at a specific node """
    
    rec_file_text = open(rec_file, "r")
    node_info = {}
    for line in rec_file_text:
        line = line.strip()
        if line.startswith("S_"):
            split_table = line.split("\t")
            node = split_table[1]
            dupes_at_node = float(split_table[2])
            trans_at_node = float(split_table[3])
            loss_at_node = float(split_table[4])
            origin_at_node = float(split_table[5])
            copies_at_node = float(split_table[6])
            node_info[node] = [os.path.basename(rec_file),dupes_at_node,trans_at_node,loss_at_node,origin_at_node,copies_at_node]
        # print(node_info)
    return node_info 

def run_node_reconstruction(uml_dir,cutoff):

    """ Function to find all the gene families at each node in the tree """

    uml_dir = Path(uml_dir)
    uml_files = list(uml_dir.glob("*.uml_rec"))
    df_final = pd.DataFrame({'Gene Family':[0], 'Duplications':[0], 'Transfers':[0], 'Loss':[0], 'Originations':[0], 'Copies':[0]})
    # print(df_final)

    for rec_file in uml_files:
        node_info = get_copies_at_node(rec_file)
        df = pd.DataFrame.from_dict(node_info, orient='index')
        df.rename(columns = {0: 'Gene Family', 1: 'Duplications', 2: 'Transfers', 3: 'Loss', 4: 'Originations', 5: 'Copies'}, inplace=True)
        df.reset_index(inplace=True)
        df_final = pd.concat([df_final, df])
        # print(df_final)

    df_final.rename(columns = {'index': 'Node'}, inplace=True)
    print(df_final)
    df_final.to_csv("Events_at_all_nodes.csv",index=False)    # this is a kind of master file with "which genes at which nodes" information to query

    df_node_cut = df_final.loc[(df_final["Duplications"] >= float(cutoff)) | (df_final["Transfers"] >= float(cutoff)) | (df_final["Loss"] >= float(cutoff)) | (df_final["Originations"] >= float(cutoff)) | (df_final["Copies"] >= float(cutoff))]
    df_node_cut.to_csv(f"Events_at_all_nodes_above_{cutoff}.csv",index=False)

    return df_final

def main():
    df = run_node_reconstruction(uml_dir, cutoff)

if __name__ == "__main__":
    """ To get all nodes: Check the internal node orders. Remember the first and last internal node orders. 
    If the number of internal nodes (branch nodes) are from 16-30, 16 and 30 should be the first and last internal node orders.

    uml_dir : path/to/dir/with/reconciliations
    cutoff : 0-1, the event should be present in this many sampled reconciliations  """

    uml_dir = sys.argv[1]
    cutoff = sys.argv[2]
    main()
```

[gene_enrichment_stats.py]()

To get genes enriched in alkaline nodes

```python
# %%
import os
import sys
import glob, csv
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact


# %%
alkaline_node_list = [79,112,123,129,133,113,80,93,97,136,138,120, "g0009", "g0011", "g0012", "g0013", "g0014", "g0015", "g0016", "g0017", "g0029", "g0030", "g0031", "g0032", "g0033", "g0034", "g0035", 'g0036']

# %%
cutoff = 0.5        # adjust for stringency

df_full = pd.read_csv('Events_at_all_nodes.csv',header=0)
df_copies = df_full[['Node','Gene Family','Copies']]
df_copies_filt = df_copies.loc[df_copies['Copies'] >= cutoff]      
df_copies_filt['alkaline_or_not'] = np.where(np.isin(df_copies_filt['Node'],alkaline_node_list), "alkaline", "not_alkaline")

# %%
df_copies_filt

# %%
# Creating a gene matrix

gene_matrix = df_copies_filt.pivot(index='Node', columns='Gene Family', values='Copies').fillna(0)
print(gene_matrix)


# %%

# convert to binary data
binary_matrix = (gene_matrix > 0).astype(int)
# print(binary_matrix)

# %%
# Creating a labels series

labels = df_copies_filt[['Node','alkaline_or_not']].drop_duplicates(['Node','alkaline_or_not']).drop('Node', axis=1)
labels = labels.squeeze()
labels.index = binary_matrix.index
print(labels)

# %% [markdown]
# ## Fisher's Exact Test
# 
# For testing significance of gene enrichment/depletion
# - Only binary data - have to turn counts to binary

# %% [markdown]
# ### Run the Fisher's test
# 
# For each gene family in a loop and then convert to dataframe
# 
# #### Interpretation of the Odds Ratio (OR)
# **OR = 1:** This indicates no association between the gene family and the genome type (alkaline vs. non-alkaline). The odds of presence in alkaline genomes are the same as the odds of presence in non-alkaline genomes.
# 
# **OR > 1:** This indicates a positive association between the gene family and alkaline genomes. The gene family is more likely to be present in alkaline genomes compared to non-alkaline genomes.
# The further the OR is greater than 1, the stronger the association. For example, an OR of 2 means that alkaline genomes are twice as likely to have the gene family present compared to non-alkaline genomes.
# 
# **OR < 1:** This indicates a negative association between the gene family and alkaline genomes, meaning the gene family is more likely to be absent in alkaline genomes compared to non-alkaline genomes.
# The further the OR is less than 1, the stronger the depletion. For example, an OR of 0.5 means that alkaline genomes are half as likely to have the gene family present compared to non-alkaline genomes.
# 
# #### Interpretation of the p-value
# p-value < 0.05 : significant
# p-value > 0.05 : non-significant
# 
# A p-value > 0.05 means that the test did not find a statistically significant difference. It does not imply that the gene family is depleted or enriched.
# 
# To determine depletion, check the odds ratio (OR). If the OR is less than 1, it suggests depletion, but the p-value must be low (typically < 0.05) for it to be statistically significant.
# 

# %%
# Initialize results
results = []

# Loop through each gene family
for gene_family in binary_matrix.columns:
    # Create contingency table
    alkaline_present = binary_matrix.loc[labels == 'alkaline', gene_family].sum()
    non_alkaline_present = binary_matrix.loc[labels == 'not_alkaline', gene_family].sum()
    alkaline_absent = len(labels[labels == 'alkaline']) - alkaline_present
    non_alkaline_absent = len(labels[labels == 'not_alkaline']) - non_alkaline_present

    contingency_table = [
        [alkaline_present, non_alkaline_present],
        [alkaline_absent, non_alkaline_absent]
    ]
    # Perform Fisher's Exact Test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
     # Store results
    results.append({
        'GeneFamily': gene_family,
        'OddsRatio': odds_ratio,
        'PValue': p_value
    })

# Convert results into a DataFrame for easy inspection
results_df = pd.DataFrame(results)

# Display the results
print(results_df)


# %% [markdown]
# ### ADJUSTING THE PVALUE USING BENJAMINI-HOCHBERG PROCEDURE for FDR (LESS CONSERVATIVE THAN BONFERRONI)
# 

# %%
# Assume `results_df` is the DataFrame containing the original p-values
# that you obtained from Fisher's Exact Test (with 'PValue' column)

# Step 1: Sort p-values in ascending order and assign ranks
results_df['Rank'] = results_df['PValue'].rank(method='min')

# Step 2: Calculate the adjusted p-values using the BH procedure formula
m = len(results_df)  # Total number of tests (gene families)
results_df['BH_Adjusted_PValue'] = results_df['PValue'] * m / results_df['Rank']

# Step 3: Ensure adjusted p-values do not exceed 1
results_df['BH_Adjusted_PValue'] = np.minimum(results_df['BH_Adjusted_PValue'], 1)

# Step 4: Determine which tests are significant after BH correction
results_df['Significant'] = results_df['BH_Adjusted_PValue'] < 0.05

# Display the results
print(results_df)

results_df_filtered = results_df.loc[results_df['PValue'] < 0.05]
print(results_df_filtered)

# %% [markdown]
# Dunno, I am not convinced by this as it is classifying some genes which I know are unique in alkaline as non-significant

# %% [markdown]
# ### Plotting the significant enriched/depleted genes in alkaline genomes

# %%
import matplotlib.pyplot as plt

# Assume 'results_df' contains the adjusted p-values and log-fold change (or odds ratio)
# If you don't have log-fold change, you can calculate the odds ratio and take the log2

# Example of adding log-fold change (using Odds Ratio or other measure of effect size)
# results_df['log2_fold_change'] = np.log2(results_df['OddsRatio'])

# Volcano plot
plt.figure(figsize=(8, 6))
plt.scatter(
    results_df['OddsRatio'],  # x-axis: log-fold change (enrichment or depletion)
    -np.log10(results_df['BH_Adjusted_PValue']),  # y-axis: -log10(p-value)
    c=results_df['BH_Adjusted_PValue'] < 0.05,  # Color points by significance (adjusted p-value < 0.05)
    cmap='coolwarm',  # Color map
    edgecolor='k',
    alpha=0.7
)

# Add significance threshold line
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='Significance Threshold')

# Labels and title
plt.title('Volcano Plot: Gene Family Enrichment/Depletion')
plt.xlabel('Odds Ratio (Enrichment/Depletion)')
plt.ylabel('-Log10 Adjusted P-value')
plt.legend()

# Show plot
plt.show()


# %% [markdown]
# ## Poisson regression
# Takes into account count data - non binary
# 

# %%
# For one gene family only

import statsmodels.api as sm

# Example: Poisson regression for count data
X = gene_matrix['pmsf_0_MSA_trimmed.ufboot.ale.uml_rec']  # Count data (number of occurrences of a gene family)
y = (labels == 'alkaline').astype(int)  # Binary outcome: Alkaline vs. Non-Alkaline


# Add constant (intercept)
X = sm.add_constant(X)


# Poisson regression model
poisson_model = sm.GLM(y, X, family=sm.families.Poisson())
poisson_results = poisson_model.fit()

print(poisson_results.summary())


# %%
import statsmodels.stats.multitest as mt

# Assuming `gene_matrix` is your dataframe of gene family counts
# and `labels` contains the corresponding 'Alkaline' or 'Non-Alkaline' labels

# Initialize a list to store p-values for each gene family
results = []

# Loop through each gene family
for gene_family in gene_matrix.columns:
    # Create the dependent variable (alkaline or not) 
    # assuming 'alkaline_or_not' is a binary column in `labels` (0 = Non-Alkaline, 1 = Alkaline)
    X = gene_matrix[gene_family]
    X = sm.add_constant(X)  # Add intercept to the model
    
    y = (labels == 'alkaline').astype(int)  # Convert labels to binary: 1 for Alkaline, 0 for Non-Alkaline
    
    # Fit Poisson regression model
    poisson_model = sm.GLM(y, X, family=sm.families.Poisson(), link=sm.families.links.log()).fit()
    
    # Get the p-value for the gene family
    results.append({
        'GeneFamily': gene_family,
        'OddsRatio': np.exp(poisson_model.params[gene_family]),
        'PValue': poisson_model.pvalues[gene_family]
    })


# %%

# Display the results
results_df_pois = pd.DataFrame(results)

# Step 1: Sort p-values in ascending order and assign ranks
results_df_pois['Rank'] = results_df_pois['PValue'].rank(method='min')

# Step 2: Calculate the adjusted p-values using the BH procedure formula
m = len(results_df_pois)  # Total number of tests (gene families)
results_df_pois['BH_Adjusted_PValue'] = results_df_pois['PValue'] * m / results_df_pois['Rank']

# Step 3: Ensure adjusted p-values do not exceed 1
results_df_pois['BH_Adjusted_PValue'] = np.minimum(results_df_pois['BH_Adjusted_PValue'], 1)

# Step 4: Determine which tests are significant after BH correction
results_df_pois['Significant'] = results_df_pois['BH_Adjusted_PValue'] < 0.05

# Display the results
print(results_df_pois)

# filtered cos it gave me some really weird odds ratios
results_df_pois_filt = results_df_pois.loc[(results_df_pois['OddsRatio'] < 50) & (results_df_pois['OddsRatio'] > 0.01)]#& ] (results_df['BH_Adjusted_PValue'] < 0.05)
print(results_df_filtered) 

# Optional: Identify significant gene families (adjusted p-value < 0.05)
# significant_gene_families = results[results['PValue'] < 0.05]
# print("Significant Gene Families:")
# print(significant_gene_families)


# %%
# Volcano plot
plt.figure(figsize=(8, 6))
plt.scatter(
    results_df_pois_filt['OddsRatio'],  # x-axis: log-fold change (enrichment or depletion)
    -np.log10(results_df_pois_filt['BH_Adjusted_PValue']),  # y-axis: -log10(p-value)
    c=results_df_pois_filt['BH_Adjusted_PValue'] < 0.05,  # Color points by significance (adjusted p-value < 0.05)
    cmap='coolwarm',  # Color map
    edgecolor='k',
    alpha=0.7
)

# Add significance threshold line
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='Significance Threshold')

# Labels and title
plt.title('Volcano Plot: Gene Family Enrichment/Depletion by Poisson Regression')
plt.xlabel('Odds Ratio (Enrichment/Depletion)')
plt.ylabel('-Log10 BH Adjusted P-value')
plt.legend()

# Show plot
plt.show()

# %%
# from Fisher's test

results_df_filt = results_df.loc[results_df['BH_Adjusted_PValue'] < 0.05]
results_df_filt

# from Poisson regression
results_df_pois_filt_2 = results_df_pois_filt.loc[results_df_pois_filt['BH_Adjusted_PValue'] < 0.05]
results_df_pois_filt_2

# %%
# merged the two tests

df_test_merged = pd.merge(results_df,results_df_pois_filt, how='outer', on='GeneFamily', suffixes=['_Fisher', '_Poisson'])
print(df_test_merged)
df_test_merged = df_test_merged.loc[(df_test_merged['BH_Adjusted_PValue_Fisher'] < 0.05) | (df_test_merged['BH_Adjusted_PValue_Poisson'] < 0.05)]
df_test_merged
# df_test_merged.drop(['Rank','BH_Adjusted_PValue','Significant'], inplace=True, axis=1)
df_test_merged
temp = df_test_merged['GeneFamily'].str.split("_", expand=True)[1]

# df_test_merged['cluster_id'] = temp['1']
df_test_merged = df_test_merged.join(temp)
df_test_merged.rename(columns={1 : 'cluster_id'},inplace=True)
df_test_merged.set_index('cluster_id', inplace=True)
df_test_merged.index = df_test_merged.index.astype(int)

df_test_merged



# %%
df_annotations = pd.read_csv('homologues.tsv', sep='\t')
df_annotations.rename(columns={'cluster id': 'cluster_id'},inplace=True)
df_annotations.set_index('cluster_id', inplace=True)
df_annotations

# %%
df_copies_filt_grouped = df_copies_filt.groupby('Gene Family').agg({'Node': lambda x: x.tolist()})
df_copies_filt_grouped

# %%
df_final = df_test_merged.join(df_annotations, how='left')
df_final

df_final.reset_index(inplace=True)
# print(df_final)
df_final.to_csv('genes_enriched_alkaline_statistics_incl_genomes.csv', index=False)

df_final = df_final[['GeneFamily','OddsRatio_Fisher','PValue_Fisher','BH_Adjusted_PValue_Fisher','OddsRatio_Poisson','PValue_Poisson','BH_Adjusted_PValue_Poisson','cluster_id','type','annotation', 'length','selective pressure (low = purifying)','codon bias',"representation", 'count','% id']]

df_final.rename({'selective pressure (low = purifying)':'selective_pressure','codon bias':'codon_bias','% id':'pid'},inplace=True)
df_final = df_final.merge(df_copies_filt_grouped, left_on='GeneFamily', right_on='Gene Family', how='left')
print(df_final)
df_final.to_csv('genes_enriched_alkaline_statistics_short.csv', index=False)
```

