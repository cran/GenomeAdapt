\name{PlotAdmix}
\alias{PlotAdmix}
\title{Plot Ancestry Admixture
}
\description{Plot admixture proportions based on their ancestries.
}
\usage{
PlotAdmix(propmat, group = NULL, col = NULL,
multiplot = TRUE, xlab = "Individuals", ylab = "Ancestry Proportion",
showgrp = TRUE, shownum = TRUE, ylim = TRUE, na.rm = TRUE)
}
\arguments{
  \item{propmat}{a sample-by-ancestry matrix of proportion estimates, returned from AdmixProp()
}
  \item{group}{a character vector of a factor according to the samples in propmat
}
  \item{col}{specify colors
}
  \item{multiplot}{single plot or multiple plots
}
  \item{xlab}{The title for the x axis
}
  \item{ylab}{The title for the y axis
}
  \item{showgrp}{show group names in the plot
}
  \item{shownum}{TRUE: show the number of each group in the figure
}
  \item{ylim}{TRUE: y-axis is limited to [0, 1]; FALSE: ylim <- range(propmat); a 2-length numeric vector: ylim used in plot()
}
  \item{na.rm}{TRUE: remove the sample(s) according to the missing value(s) in group
}
}
\details{The minor allele frequency and missing rate for each SNP passed in snp.id are calculated over all the samples in sample.id.
}
\value{ Return to an anscrstry proportion plot/population structure plot that the same as STRUCTURE plot.
}
\references{
Zheng X, Weir BS. Eigenanalysis on SNP Data with an Interpretation of Identity by Descent. Theoretical Population Biology. 2015 Oct 23. pii: S0040-5809(15)00089-1. doi: 10.1016/j.tpb.2015.09.004.

Zheng, X., Levine, D., Shen, J., Gogarten, S. M., Laurie, C., & Weir, B. S. (2012). A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics, 28(24), 3326-3328.
}

\examples{
##---- Do genome scan and get the ancestry proportation for individuals ----
# Scan genomes (HapMap)

HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),
method="EIGMIX",num.thread = 1L, autosome.only=TRUE,
remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)


# get population information
genofile <- SNPRelate::snpgdsOpen(SNPRelate::snpgdsExampleFileName())
pop_code <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
SNPRelate::snpgdsClose(genofile)

# define groups
groups <- list(CEU = samp.id[pop_code == "CEU"],
    YRI = samp.id[pop_code == "YRI"],
    CHB = samp.id[is.element(pop_code, c("HCB", "JPT"))])

### estimate the ancestry proportion

Admixpro=AdmixProp(HapmapScan,groups=groups,bound=TRUE)

PlotAdmix(Admixpro,group=as.factor(pop_code),multiplot = FALSE)
}
\keyword{plot}

