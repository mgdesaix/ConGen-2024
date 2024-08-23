# ConGen-2024: Population assignment with genotype likelihoods 

This is a walk-through of performing population assignment with genotype likelihoods from low-coverage whole genome sequencing (lcWGS) using the software [WGSassign](https://github.com/mgdesaix/wgsassign). The materials here correspond to the activities in my hands-on session at [ConGen 2024](https://www.umt.edu/ces/conferences/congen/), and ideally will doubly serve as an augmented vignette of the instructions for WGSassign. 

## Pre-reqs

These activities assume the reader has an understanding of R (and tidyverse!) and command line coding, as well as access to a Linux system for running WGSassign (Note: If you are using the ConGen RStudio server, it already has WGSassign installed to a Conda environment (`conda activate WGSassign`)! Thanks IT!). Many of the files are also provided in the repository to run through the R code.

## Overview of activities

This session will be broken into the following components:

1.  [Working with genotype likelihood data: 01-genotype-likelihood-data.md](./01-genotype-likelihood-data.md)

2.  [WGSassign basics: 02-wgsassign-basics.md](./02-wgsassign-basics.md)

3.


## Related resources

-  WGSassign software: <https://github.com/mgdesaix/wgsassign>
-  Initial WGSassign publication: [Population assignment from genotype likelihoods for low-coverage whole-genome sequencing data](https://doi.org/10.1111/2041-210X.14286)
-  Example of use: [Low-coverage whole genome sequencing for highly accurate population assignment: Mapping migratory connectivity in the American Redstart (Setophaga ruticilla)](https://doi.org/10.1111/mec.17137)
