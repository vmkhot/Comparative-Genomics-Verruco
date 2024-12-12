
## Orthologs
Homologous proteins were clustered using mmseq2 
 >Homologous proteins and RNA genes are clustered using [mmseqs](https://github.com/soedinglab/MMseqs2) version 15-6f452 with --min-seq-id 0.5

Orthologs and paralogs are called
>Orthologues and paralogues are called based on the distance to the center-protein

Was a total of ~26926 groups from CDS

296 accepted - RNA

Filtering:
Groups with less 3 representative taxa were removed + 100% paralogs - 12681 removed
    - 14245 CDS groups remained
Orthologous groups with less than 4 sequences are not used as iqtree cannot produce bootstraps for them
 - this removed ~6360 out of 14245
 - **7885 groups remained**
- *These were used for the multiple sequence alignment*


## MSA

```
[vmkhot@fc106 30_CDS]$ ls -1 seq_mt_75/ | wc -l 
1
[vmkhot@fc106 30_CDS]$ ls -1 seq_lt_25 | wc -l 
102
[vmkhot@fc106 30_CDS]$ ls -1 seq_50_to_75 | wc -l 
9
[vmkhot@fc106 30_CDS]$ ls -1 seq_25_to_50/ | wc -l 
18
```

Famsa is used for the multiple sequence alignment

> Each cluster of homologous proteins is aligned with [famsa](https://github.com/refresh-bio/FAMSA) version 2.2.2

## Trimming MSA

trimAl is used for trimming the MSA using the parameters `-gappyout` 

## Gene Trees

Gene trees were made using iqtree2 using the pmsf tree inference and the model LG+C20+G. 192 trees were remade by using the option -keep-ident so that non-unique sequences aren't discarded from the analysis, leaving < 4 sequences behind. This likely affects the bootstrapping but don't know exactly how.

```bash
iqtree2 -s $f -m LG+F+G -T 5 --mem 2G -pre seq_lt_25/dir_103/guide_trees/pmsf_guide_${newname}

iqtree2 -s $f -m LG+C20+G -ft seq_lt_25/dir_103/guide_trees/pmsf_guide_${newname}.treefile -keep-ident -mwopt -B 5000 -wbtl -T 5 --mem 2G -pre seq_lt_25/dir_103/pmsf_${newname}
```
