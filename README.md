# CircSeqAlignTk

CircSeqAlignTk is designed for end-to-end RNA-Seq data analysis
of circular genome sequences, from alignment to visualization.
It mainly targets viroids which are composed of 246-401 nt circular RNAs.
In addition, CircSeqAlignTk implements a tidy interface
to generate synthetic sequencing data that mimic real RNA-Seq data,
allowing developers to evaluate the performance of alignment tools
and workflows.


## Installation

To install the CircSeqAlignTk package,
start R (≥ 4.2) and run the following steps:

```
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install("CircSeqAlignTk")
```


## Documentation

Vignette, including overview of CircSeqAlignTk package,
detailed explanation of functions, and several case studies,
can be accesssed with the following steps:

```
browseVignettes('CircSeqAlignTk')
```


