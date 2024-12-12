**GTDB CLASSIFICATION**: (FAMILY LEVEL) d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Opitutales;f__T3Sed10-336;g__;s__
**Plan**
- soda lake opitutales (genome of interest) ( + 1)
- 2 other soda lake opitutales (from JZ paper) (+ 2)
    - GCA_007692545.1_ASM769254v1
    - GCA_007695295.1_ASM769529v1
- All genomes in the family T3Sed10-336 (+ 15)
- 1 family-lvl representative from the rest of the families? (+ 27)
- 1 order-lvl rep from remainder of orders? (+ 7)
- 1 class-lvl rep "" (+ 3)
- 5 random genomes from Plancto phylum as outgroup (+ 5)

60 genomes in total
* From the gtdb metadata file, there's a column for "isolation source" - to check where the genomes were isolated from

**In reality - there are 65 genomes in total:**
- soda lake opitutales (genome of interest) ( + 1)
- 2 other soda lake opitutales (from JZ paper) (+ 2)
- All genomes in the family T3Sed10-336 (+ 10)
    - these were less because of completeness filters
- 2 family-lvl representative from the rest of the families? (+ 37)
    - doubled number of family level reps
- 1 order-lvl rep from remainder of orders? (+ 7)
- 1 class-lvl rep "" (+ 3)
- 5 random genomes from Plancto phylum as outgroup (+ 5)


After running metaerg - I found 7 genomes that had to be discarded because they are not complete. < 1000 genes were annotated on each. Need to use retro-check completeness (instead of checkm2 for this)

Decided to redo the genome selection to include more representatives, because 56 genomes in total did not feel like enough representation from within the Verrucomicrobiaceae order

- soda lake opitutales (genome of interest) ( + 1)
- 2 other soda lake opitutales (from JZ paper) (+ 2)
- All genomes in the family T3Sed10-336 (+ 10)
    - these were less because of completeness filters
- 3 family-lvl representative from the rest of the families? (+ 39)
    - doubled number of family level reps
- 2 order-lvl rep from remainder of orders? (+ 14)
- 1 class-lvl rep "" (+ 3)
- 5 random genomes from Plancto phylum as outgroup (+ 5)

