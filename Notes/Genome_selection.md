In the case of Opitutales,

**GTDB CLASSIFICATION**: (FAMILY LEVEL) d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Opitutales;f__T3Sed10-336;g__;s__

**CLOSEST GTDB PLACEMENT** d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Opitutales;f__T3Sed10-336;g__T3Sed10-336;s__T3Sed10-336 sp003566375, GCA_003566375.1 at 77% ANI. 

**FILTERED NCBI TAXONOMY** d__Bacteria; p__Verrucomicrobia; c__Opitutae; o__; f__; g__; s__ 

**UNFILTERED NCBI TAXONOMY** d__Bacteria; x__PVC group; p__Verrucomicrobia; c__Opitutae; x__unclassified Opitutae; s__Opitutae bacterium

OPTIONS

1. Everything from GTDB family level and lower from f__T3Sed10-336 --> 12 genomes --> not enough resolution for genome species tree? 
2. Everything from filtered NCBI taxonomy class level and lower of c__Opitutae --> 1400+ genomes -->is this too many?
3. Everything from unfiltered NCBI taxonomy species level and lower of s__Opitutae bacterium --> 59 genomes
4. All soda lake and lab bioreactor MAGs associated with Opitutales genomes.

Tried option #2 MCL quit, out of memory handler at 100G. Opened the gtdbtk tree - where to chop this? Is it okay to use "representative" genomes from this tree? - enough resolution?


Pick the whole family from GTDB tree and then about 15 from the sister clade (GCFs)

Picked 44 other genomes to compare with the Opitutales genome, including  all of them from the same clade (blue) + refseq ones from the sister clade + extra random genbank genomes from sister clade. 

1 genome came from the outgroup to root the tree later - include a couple more for the outgroup

[[genomes_selected_tree.png]]

## CheckM2

Have to test whether the genomes are complete using CheckM2
 #completeness

```bash
checkm2 predict -t 30 -x fna --input ./genomes_1 --output-directory ./Checkm2_1
```

[[checkm2 results.xlsx]]  
Sheet 2 has checkm2 results from all the genomes. Only genomes > 90% completeness were selected for further analyses.

## FastANI

#ani

FastANI is to assess average nucleotide identity between genomes. (How similar are they on a whole genome level?)

```bash
source ~/miniconda3/etc/profile.d/conda.sh

conda activate fastani

fastANI --ql all_genomes.list --rl all_genomes.list -o fastani_results.tsv --fragLen 1000
```

Results  
Only similarities >80% are significant.

```bash
# to filter 
awk '{if ($3 >= 80) print $0}' fastani_results.tsv > fastani_results_filtered.tsv
```
[[fastani_results_filtered.txt]]

Finally, 47 genomes were processed  
Next, we have to run [MetaErg](https://github.com/kinestetika/MetaErg) to produce annotated amino acid files. #metaerg