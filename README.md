# scRNAseq_rawls5532
Collection of work for the single-cell RNA sequencing dataset described in the following: https://doi.org/10.1101/2020.12.13.422569

This will include two functions (that will be periodically updated) that were used as part of custom scripts that interface with the packages `tidyverse` and `Seurat`. I will also include code used to process the data, and some examples of figures generated.


The main file is the `WT Reference Pickout... .R` and it has been roughly partitioned into several sections. For the most part it follows the `Seurat` vignette on reference sample-based integration of single-cell data but it deviates in subsequent analyses. I show here examples of "painting" cells of a certain clustering from one UMAP to another. At the time when I did this, there was not really a straightforward way to do this although some additional functionalities have been added to `Seurat` that make this easier to do out of the box. I also include analyses that may not have made it into the paper either due to publication size limits or perhaps they were meant to be more exploratory in nature and didn't really answer our questions directly.
