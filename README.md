# CircSeqAlignTk

*CircSeqAlignTk* is an end-to-end toolkit for RNA-Seq data analysis of
circular genome sequences, supporting tasks from alignment to visualization.
It is primarily designed for studying viroids, which are small, circular RNAs
typically composed of a few hundred nucleotides.
In addition to analysis capabilities, *CircSeqAlignTk* provides a streamlined
interface for generating synthetic RNA-Seq data that closely mimics real datasets.
This functionality enables developers to benchmark alignment tools,
test novel alignment algorithms, and validate new workflows.
This vignette offers an overview of *CircSeqAlignTk* and demonstrates its usage.


# Installation

To install the *CircSeqAlignTk* package, start R (≥ 4.2) and run the following steps:

```{r install_package, eval=FALSE}
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('CircSeqAlignTk')
```

Note that to install the latest version of the *CircSeqAlignTk* package,
the latest version of R is required.



## Documentation

Vignette, including overview of *CircSeqAlignTk* package,
detailed explanation of functions, and several case studies,
can be accessed with the following code.

```
browseVignettes('CircSeqAlignTk')
```


## Citation

Sun J, Fu X and Cao W.
CircSeqAlignTk: An R package for end-to-end analysis of RNA-seq data for circular genomes.
F1000Research 2024, 11:1221.
doi: [10.12688/f1000research.127348.2](https://doi.org/10.12688/f1000research.127348.2).



