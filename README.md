[![DOI](https://zenodo.org/badge/409003592.svg)](https://zenodo.org/badge/latestdoi/409003592)

# scTfh
This repository contains the code used in our single cell sequencing paper: "Follicular T Cells are Clonally and Transcriptionally Distinct in B Cell-Driven Autoimmune Disease"

## Installation guide
No installation necessary

## Demo
Model data is included in the `data` directory.

## Instructions
1. Download `data` and `code` directories
2. Set `data` as the working directory
3. Download `gex` and `vdj` data from GSE157649, available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157649
4. Run each script in numerical order in either Rstudio or Python.  These scripts will generate the figures presented in our manuscript.\
\
NB: Expected total run time is 3-5 days

## System requirements and software
CellRanger	10X Genomics	v4.0.0\
R	The Comprehensive R Archive Network	v3.6.1\
Rstudio v1.2.1578\
2017PWM	https://bitbucket.org/slofgren/antigen-id-paper-code/src/master/	N/A\
GLIPH2	http://50.255.35.37:8080/	v2\
scanpy	https://pypi.org/project/scanpy/	v1.5.1\
vegan	https://cran.r-project.org/web/packages/vegan/	v2.5-6\
qgraph	https://cran.r-project.org/web/packages/qgraph	v1.6.5\
immunarch	https://immunarch.com/	v0.5.5\
ggseqlogo	https://omarwagih.github.io/ggseqlogo/	v0.1\
msigdbr	https://cran.r-project.org/web/packages/msigdbr	v7.0.1\
clusterProfiler	https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html	v3.14.3\
pheatmap	https://cran.r-project.org/web/packages/pheatmap/	v1.0.12\
vennDiagram	https://cran.r-project.org/web/packages/VennDiagram	v1.6.20\
ggpubr	https://rpkgs.datanovia.com/ggpubr/	v0.2.5.999\
EnhancedVolcano	https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html	v1.4.0\
biomaRt	https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html	v2.42.0\
ggplot2	https://ggplot2.tidyverse.org/	v3.3.0\
Seurat	https://satijalab.org/seurat/	v3.1.4\
anndata2ri	https://github.com/theislab/anndata2ri	v1.0.2\
DEseq2	https://bioconductor.org/packages/release/bioc/html/DESeq2.html	v1.26.0\
DeepTCR	https://github.com/sidhomj/DeepTCR	v1.4.15


