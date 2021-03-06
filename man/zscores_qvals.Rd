\name{zscores_qvals}
\alias{zscores_qvals}
\title{Calculating (q)p-values from multiple factors/scores
}
\description{ Converting the Z-scores to (q)p-values. This function calibrates the p-value/q-values considering multiple scores based on Mahalanobis Distance/componentwise method.
}
\usage{
zscores_qvals(x, outlier.method = "mahalanobis",
estim = "pairwiseGK", pval.coret.method = "bonferroni")
}
\arguments{
  \item{x}{ The GenomeAdapt object containing the locus scores}
  \item{outlier.method}{ The methods used to detect the outliers, "mahalanobis" or "componentwise", default is mahalanobis}
  \item{estim}{ Method used to estimate Mahalanobis distance. The choices are : "mcd" for the Fast MCD algorithm of Rousseeuw and Van Driessen, "weighted" for the Reweighted MCD, "donostah" for the Donoho-Stahel projection based estimator, "M" for the constrained M estimator provided by Rocke, "pairwiseQC" for the orthogonalized quadrant correlation pairwise estimator, and "pairwiseGK" for the Orthogonalized Gnanadesikan-Kettenring pairwise estimator. The default "auto" selects from "donostah", "mcd", and "pairwiseQC" with the goal of producing a good estimate in a reasonable amount of time.
}
  \item{pval.coret.method}{ Correction method for p-values,choices are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none".The adjustment methods include the Bonferroni correction ("bonferroni") in which the p-values are multiplied by the number of comparisons. Less conservative corrections are also included by Holm (1979) ("holm"), Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"), Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), and Benjamini & Yekutieli (2001) ("BY"), respectively. A pass-through option ("none") is also included.

}
}
\details{Calculating (q)p-values from GenomeAdapt to identify the outlier loci

}
\value{ A dataframe with p-values, adjusted p-values and q-values.

\item{pvals}{ A dataframe containing 3 columns, including the p-values, q-values, and adjusted p-values for all loci }
\item{chr }{The chromsomes for the dataset}

}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289-300. http://www.jstor.org/stable/2346101.

R. A. Maronna and R. H. Zamar (2002) Robust estimates of location and dispersion of high-dimensional datasets. Technometrics 44 (4), 307-317.

Capblancq, T., Luu, K., Blum, M. G., & Bazin, E. (2018). Evaluation of redundancy analysis to identify signatures of local adaptation. Molecular Ecology Resources, 18(6), 1223-1233.

}
\examples{
##---- Do genome scan ----
HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),
method="EIGMIX",num.thread = 1L, autosome.only=TRUE,
remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

## estimating the q-values from genome scan

\donttest{Hapmapqval=zscores_qvals(HapmapScan)}

}

\keyword{ p-values }
