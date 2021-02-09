# Loop-analysis-scripts
## Full analysis pipeline
NOTE: Keep all files (including .bedpe files and hg38 gff file) in a single directory. For the python scripts and Linux terminal, use this as the working directory.

Dependencies: bedtools, Python >= 3.7, pandas >= 1.2.0 

To replicate this analysis,

1) Run the Common_Loop_generator.py script (preferably in your python IDE). The output should be a set of .bed files.
2) Run the following in your Linux terminal:
```
bedtools window -a anyLoops_AML_1.bed -b hg38.gff -w 15000 > mapped_anyLoops_AML_1.csv
bedtools window -a anyLoops_AML_2.bed -b hg38.gff -w 15000 > mapped_anyLoops_AML_2.csv
bedtools window -a anyLoops_Knee_1.bed -b hg38.gff -w 15000 > mapped_anyLoops_Knee_1.csv
bedtools window -a anyLoops_Knee_2.bed -b hg38.gff -w 15000 > mapped_anyLoops_Knee_2.csv
bedtools window -a blacklistLoops_1.bed -b hg38.gff -w 15000 > mapped_blacklistLoops_1.csv
bedtools window -a blacklistLoops_2.bed -b hg38.gff -w 15000 > mapped_blacklistLoops_2.csv
bedtools window -a commonLoops_AML_1.bed -b hg38.gff -w 15000 > mapped_commonLoops_AML_1.csv
bedtools window -a commonLoops_AML_2.bed -b hg38.gff -w 15000 > mapped_commonLoops_AML_2.csv
bedtools window -a commonLoops_Knee_1.bed -b hg38.gff -w 15000 > mapped_commonLoops_Knee_1.csv
bedtools window -a commonLoops_Knee_2.bed -b hg38.gff -w 15000 > mapped_commonLoops_Knee_2.csv
bedtools window -a middleLoops_AML_1.bed -b hg38.gff -w 15000 > mapped_middleLoops_AML_1.csv
bedtools window -a middleLoops_AML_2.bed -b hg38.gff -w 15000 > mapped_middleLoops_AML_2.csv
bedtools window -a middleLoops_Knee_1.bed -b hg38.gff -w 15000 > mapped_middleLoops_Knee_1.csv
bedtools window -a middleLoops_Knee_2.bed -b hg38.gff -w 15000 > mapped_middleLoops_Knee_2.csv
```
3) Run the Get_genes.py script (preferably in your python IDE). The output should be a set of .txt query files containing genes associated with the following strata of loops:
```
query_genes_mapped_anyLoops_AML.txt #loops specific to AML found in any sample
query_genes_mapped_anyLoops_Knee.txt #loops specific to Knee found in any sample
query_genes_mapped_commonLoops_AML.txt #loops specific to AML found in all 3 samples
query_genes_mapped_commonLoops_Knee.txt #loops specific to Knee found in all 3 samples
query_genes_mapped_blacklistLoops.txt #loops found shared between AML and Knee

query_genes_mapped_middleLoops_AML.txt #anyLoops_AML, complement commonLoops_AML
query_genes_mapped_middleLoops_Knee.txt #anyLoops_Knee, complement commonLoops_Knee
```

## Query file output format

| gene | name |	biotype |	description |	feature_chr |	feature_start |	feature_end |	strand |	Name |	Role in Cancer |	ind |	fdrBL |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| ENSG00000116731	| PRDM2	| protein_coding	| PR/SET domain 2	| 1	| 13700198.0 | 13825079.0 |	+	| PR/SET domain 2| TSG | AML29_9431; Knee50_6156 |	0.00016261874; 2.3308723e-07|
| ENSG00000132906	| CASP9	| protein_coding	| caspase 9	| 1	| 15490832.0	| 15526534.0	| -	| caspase 9	| TSG	| AML28_12698; Knee49_5917	| 0.0011177248; 0.004380403
| ENSG00000117118	| SDHB	| protein_coding	| succinate dehydrogenase complex iron sulfur subunit B	| 1	| 17018722.0	| 17054032.0	| -	| succinate dehydrogenase complex, subunit B, iron sulfur (Ip)	| TSG	| AML28_12697; Knee50_6825; AML28_11210; AML30_4744; Knee49_6057; Knee50_6021	| 0.0011177248; 0.0031032744; 7.513726e-06; 0.003536516; 0.0015313101; 7.0954746e-05
| ENSG00000074964	| ARHGEF10L	| protein_coding	| Rho guanine nucleotide exchange factor 10 like	| 1	| 17539698.0	| 17697874.0	| +	| Rho guanine nucleotide exchange factor 10 like	| TSG	| AML28_11667; AML29_9294; AML30_5320; Knee47_2373; Knee49_6298; Knee50_6323	| 1.6974686e-37; 3.1449943e-34; 1.0405936e-18; 0.005654519; 3.7536597999999996e-19; 6.3534770000000005e-28
