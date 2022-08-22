# CircSeqAlignTk

CircSeqAlignTk is an R package for end-to-end analysis,
from alignment to visualization,
of RNA-seq data obtained from organelles or organisms
with circular genome sequences,
such as bacteria, viruses, and viroids.
In addition, CircSeqAlignTk implements a tidy interface
to generate synthetic sequencing data that mimic real RNA-seq data,
allowing developers to evaluate the performance of alignment tools,
new alignment algorithms, and new workflows.


## Installation

To install the CircSeqAlignTk package,
start R (≥ 4.2) and run the following steps:

```
# install the dependency packages
install.packages(c('tidyverse', 'devtools'))
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c('BiocGenerics', 'Biostrings', 'IRanges',
                       'Rbowtie2', 'Rhisat2', 'Rsamtools', 'ShortRead'))

library(devtools)
install_github('jsun/CircSeqAlignTk', build_vignettes = TRUE)
```


## Documentation

Vignette, including overview of CircSeqAlignTk package,
detailed explanation of functions, and several case studies,
can be accesssed with the following steps:

```
browseVignettes('CircSeqAlignTk')
```


