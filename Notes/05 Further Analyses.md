alkaline nodes:

94 (except g990067)
92(except g990067)

At the outermost
 112, 79, 113, 80, 129, 93, 97, 136
120 (only g0017)
81 (except g99065)
internal nodes:

123, 133, 138 (except g0024, g990075, g0054)
94 (except g99067)


|       |     |
| ----- | --- |
| g0009 |     |
| g0010 |     |
| g0011 |     |
| g0012 |     |
| g0013 |     |
| g0014 |     |
| g0015 |     |
| g0016 |     |
| g0017 |     |
|       |     |
| g0029 |     |
| g0030 |     |
| g0031 |     |
| g0032 |     |
| g0033 |     |
| g0034 |     |
| g0035 |     |
| g0036 |     |


## Using PastML for character mapping 

Used PastML for character mapping with the MPPA tree joining method and the F81 model (as is the defaults of the program). This was to identify nodes of interest. 


![[Pasted image 20241210142241.png]]

Image from PastML:
- yellow - marine
- lavender - gut
- dark blue - alkaline
- green - freshwater
- pink - soil
- orange - hot spring

Used 79,112,123,129,133,113,80,93,97,136,138,120, "g0009", "g0011", "g0012", "g0013", "g0014", "g0015", "g0016", "g0017", "g0029", "g0030", "g0031", "g0032", "g0033", "g0034", "g0035", 'g0036' are "alkaline" nodes/leaves

## Test statistics for finding genes enriched in alkaline genomes

For binary data - can use fisher's exact test or chi-squared test. Both require contingency tables (below)

Fisher's Exact Test uses a **2x2 contingency table** to compare categorical data. For each gene family, the table looks like this:

||**Alkaline**|**Non-Alkaline**|**Total**|
|---|---|---|---|
|**Gene Present**|A|B|A + B|
|**Gene Absent**|C|D|C + D|
|**Total**|A + C|B + D|N|

Where:

- `A`: Number of alkaline genomes with the gene family present.
- `B`: Number of non-alkaline genomes with the gene family present.
- `C`: Number of alkaline genomes with the gene family absent.
- `D`: Number of non-alkaline genomes with the gene family absent.

To account for raw count data, have to use Poisson regression or negative binomial (if variance > mean)

First, I'll try the fisher's test

Correcting for pvalues for mulitple testing:
1. Bonferroni
    - Multiplies the pvalue with the number of tests so for a large number of gene families - this becomes very conservative
2. The chat suggest Benjamini-Hochberg (for correcting false discovery rate)


### Results

75 gene families are enriched according to Fisher's test + BH correction in alkaline genomes

0 gene families according to Poisson + BH

the results are in the selected_genomes workbook, genes_enriched_alkaline_genomes worksheet. 

The higher the Odd's Ratio (>1) - the more enriched they are 
The lower the Odd's ratio (< 1) - the more depleted they are

Can run the jupyter notebook (gene_enrichment_stats.ipynb) again with different nodes if needed. The raw results are in the csv files that are generated:
- genes_enriched_alkaline_statistics_incl_genomes.csv
    - with genomes/genes
- genes_enriched_alkaline_statistics_short.csv
    - shorter file/summary


**Nodes that might be interesting for transitions**

There are three transitions to soda lakes from marine - one larger one at Node 142, then again at Node 94 and Node 81
There's an interesting transition from marine to gut?? - 131 and 135
