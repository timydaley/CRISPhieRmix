# CRISPhieRmix

A hierarchical mixture model for large-scale CRISPRi and CRISPRa pooled screens.  CRISPhieRmix uses a two groups model ([Efron 2008](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.334.4762&rep=rep1&type=pdf)) to compute the local false discovery, the posterior probability that a gene is null and has no effect in the screen.  CRISPhieRmix can be used with or without negative control guides, though we find that effective control of the false discovery rate can only be acheived with negative control guides due to long tails in the null distribution.  This can only be estimated with negative control guides.  
CRISPhieRmix is written in R and C++.  To use CRISPhieRmix in R, install the devtools package (https://cran.r-project.org/web/packages/devtools/index.html) and then type
```
devtools::install_github("timydaley/CRISPhieRmix")
```
CRISPhieRmix takes as input normalized log fold changes of the individual gene targeting guides, along with the associated genes (in the same order as the log fold changes), with the option of including the log fold changes of the negative control guides.  We typically use log2 fold changes computed by DESeq2, although in theory any properly normalized log fold changes can be used.  To use CRISPhieRmix in R, type
```
CRISPhieRmix::CRISPhieRmix(x = log2fc, geneIds = geneIds, negCtrl = log2fc.negCtrl)
```

Full details on how to use CRISPhieRmix, with detailed instructions and examples can be found in the manual and vignette, in the vignette section. 

If you have any questions or comments, please email me at tdaley@stanford.edu.  For any issues with the software, please create an issue through github with the Issues link above.  

CRISPhieRmix is distributed with a GPL license, for full details see the file gpl-3.0.txt.
