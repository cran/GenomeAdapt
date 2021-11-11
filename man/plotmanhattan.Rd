\name{plotmanhattan}
\alias{plotmanhattan}
\title{Making manhattan plot
}
\description{Ploting a Manhattan plot showing the (q)p-values for each SNP.
}
\usage{
plotmanhattan(x, ylim=c(0,200),xlab="",
ylab="-log(p-value)",col = x$chr, pch="*",h=10, lcol="blue", ...)
}
\arguments{
  \item{x}{The zscores_qvals object
}
  \item{ylim}{ The y axis range
}
  \item{xlab}{The x title
}
  \item{ylab}{The y title
}
  \item{col}{The colour used for indicating chromosomes
}
  \item{pch}{The shape of of the data points.
}
  \item{h}{ The cutoff line
}
  \item{lcol}{The colour of the cutoff line
}
  \item{...}{passing to other parameters of plot
}
}

\value{ Return to a manhattan plot with (q)p-values for each SNP.

}

\examples{
### using Hapmap data
HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),
method="EIGMIX",num.thread = 1L, autosome.only=TRUE,
remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

### Not running, it takes a while to finish this
\donttest{ Hapmapqval=zscores_qvals(HapmapScan)

## plot
plotmanhattan(Hapmapqval$pvals$pvals$p.values,col=Hapmapqval$chr)
}

}
\keyword{plot}

