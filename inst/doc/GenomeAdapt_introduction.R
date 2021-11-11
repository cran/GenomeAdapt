## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- library-----------------------------------------------------------------
#Install from CRAN
#install.packages("GenomeAdapt")
## or you can get the latest version of HierDpart from github
#library(devtools)
#install_github("xinghuq/GenomeAdapt")

library("GenomeAdapt")



## -----------------------------------------------------------------------------
# example genepop file
# using Hapmap data
#Do genome scan
HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),method="EIGMIX",num.thread = 6, autosome.only=TRUE, remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

## it takes a while to finish this



## ----fig1, fig.height = 5, fig.width = 8.5, fig.align = "center"--------------
# get population information
library(SNPRelate)
genofile <- SNPRelate::snpgdsOpen(SNPRelate::snpgdsExampleFileName())
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snpgdsClose(genofile)

# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

### estimate the ancestry proportion

Admixpro=AdmixProp(HapmapScan,groups=groups,bound=TRUE)

PlotAdmix(Admixpro,group=as.factor(pop_code), xlab="Individuals",multiplot = FALSE)


## ----fig2, fig.height = 10, fig.width = 10, fig.align = "center"--------------

pairs(HapmapScan$eig$vectors[,1:10], col = as.factor(pop_code), pch = 19)



## ----fig3, fig.height = 10, fig.width = 9, fig.align = "center"---------------


## Calculating the Z-score for loci  
 Hapmapqval=zscores_qvals(HapmapScan)
## plot

 plotmanhattan(Hapmapqval$pvals$pvals$p.values,col=Hapmapqval$chr)


