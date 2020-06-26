# imputation-pipeline

This repository contains scripts used in an imputation pipeline for low-coverage (0.05-1x) ancient genomes.

### imputation_pipeline.snakefile
Default pipeline, including down-sampling, genotype likelihood calling in [ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home), and imputation in [Beagle](http://faculty.washington.edu/browning/beagle/beagle.html). It is written as a [snakemake](https://snakemake.readthedocs.io) workflow, but bash commands can be easily extracted from the "shell" sections.   

### local_imputation_pipeline.sh
(Author: [Eugenia D'Atanasio](https://www.ibpm.cnr.it/index.php?option=com_cnr&view=profile&id=613&lang=en#))

Pipeline used for local imputation when only genotypes at certain positions are of interest. 

### pmd_filter.py
Filters out GLs that could result from deamination (C -> T changes).

### compare_accuracy.py
Compares imputed and true genotypes to produce the accuracy table for various allele frequency bins. 
