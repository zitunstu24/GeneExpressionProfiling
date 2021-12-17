# GeneExpressionProfiling

1. Landscape of gene expression of Tomato. 
2. Genotype specific heat responsive genes patterns 
3. tissue specific heat responsive genes patterns 
4. Candidate gene selection
5. Co-expression of candidate genes for their co regulation and Conserved candidate genes selection by phylogenetic profiling
6. Gene ontology enrichment analysis of candidate genes


# Quick start 

Following the steps

For mac users, it is recommended to install xquartz (https://github.com/XQuartz/XQuartz/releases/download/XQuartz-2.8.1/XQuartz-2.8.1.dmg) for handling libraries and applications.  

1. $ git clone https://github.com/zitunstu24/GeneExpressionProfiling
2. $ cd GeneExpressionProfiling
3. $ Rscript packages.R
4. $ cd inputs
5. $ gunzip combined_abundance.tsv.gz
6. $ cd ../

We are ready to go

Then run scripts in terminal with Rscript (e.g Rscript 1_DEGs.R) 

# Features

RNA samples collected from GEO with a diverse combination of 
1. 17 genotypes of Tomato
2. 7 tissue types
3. 3 level of heat treatments with respective controls
4. Different duration of heats

inputs/master.table.csv contains all details information of samples, SRA, genotype, tissues and stress. 

Downloading, quality checking and mapping was done by using https://github.com/NAMlab/rnaseq-mapper tool.

Co-evolution profil was analysed by following Hogprof (https://github.com/DessimozLab/HogProf) and OMA browser. 

If you have any questions please feel free to contact with me zitunstu24@gmail.com (Abul Khayer)
