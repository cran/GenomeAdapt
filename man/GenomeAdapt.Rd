\name{GenomeAdapt}
\alias{GenomeAdapt}
\alias{GenomeAdapt.bed}
\alias{GenomeAdapt.vcf}
\alias{GenomeAdapt.gds}
\title{Detecting signatures of local adaptation based on multidimensional ancestry genetic map
}
\description{This function implements genome scan to identify the signatures of local adaptation. See details.
}
\usage{
GenomeAdapt(genfile, method = "EIGMIX", sample.id = NULL, snp.id =
                    NULL, autosome.only = TRUE, remove.monosnp = TRUE, maf
                    = NaN, missing.rate = NaN, num.thread = 1L, out.fn =
                    NULL, out.prec = c("double", "single"), out.compress =
                    "LZMA_RA", with.id = TRUE, verbose = TRUE,...)

## S3 method for class 'GenomeAdapt.bed'
GenomeAdapt.bed(genfile, method = "EIGMIX", sample.id = NULL, snp.id =
                    NULL, autosome.only = TRUE, remove.monosnp = TRUE, maf
                    = NaN, missing.rate = NaN, num.thread = 1L, out.fn =
                    NULL, out.prec = c("double", "single"), out.compress =
                    "LZMA_RA", with.id = TRUE, verbose = TRUE, ...)

## S3 method for class 'GenomeAdapt.vcf'
GenomeAdapt.vcf(genfile, method = "EIGMIX", sample.id = NULL, snp.id =
                    NULL, autosome.only = TRUE, remove.monosnp = TRUE, maf
                    = NaN, missing.rate = NaN, num.thread = 1L, out.fn =
                    NULL, out.prec = c("double", "single"), out.compress =
                    "LZMA_RA", with.id = TRUE, verbose = TRUE, ...)

## S3 method for class 'GenomeAdapt.gds'
GenomeAdapt.gds(genfile, method = "EIGMIX", sample.id = NULL, snp.id =
                    NULL, autosome.only = TRUE, remove.monosnp = TRUE, maf
                    = NaN, missing.rate = NaN, num.thread = 1L, out.fn =
                    NULL, out.prec = c("double", "single"), out.compress =
                    "LZMA_RA", with.id = TRUE, verbose = TRUE, ...)
}

\arguments{
  \item{genfile}{ Genotype file containg sample ID and SNP ID. Genotype format can be plink, vcf, or GDS file.

}
\item{method}{The method used to measure IBD. Default is "EIGMIX" according to Zheng, X., & Weir, B. S. (2016). "GCTA" - genetic relationship matrix defined in CGTA; "Eigenstrat" - genetic covariance matrix in EIGENSTRAT; "EIGMIX" - two times coancestry matrix defined in Zheng & Weir (2015), "Weighted" - weighted GCTA, as the same as "EIGMIX", "Corr" - Scaled GCTA GRM (dividing each i,j element by the product of the square root of the i,i and j,j elements), "IndivBeta" - two times individual beta estimate relative to the minimum of beta.

}
\item{sample.id}{a vector of sample id specifying selected samples; if NULL, all samples are used

}
\item{snp.id}{a vector of snp id specifying selected SNPs; if NULL, all SNPs are used

}
\item{autosome.only}{use autosomal SNPs only; if it is a numeric or character value, keep SNPs according to the specified chromosome

}
\item{remove.monosnp}{remove monomorphic SNPs

}
\item{maf}{filter SNPs with ">= maf" only; if NaN, no MAF threshold

}
\item{missing.rate}{filter the SNPs with "<= missing.rate" only; if NaN, no missing threshold

}
\item{num.thread}{the number of (CPU) cores used; if NA, detect the number of cores automatically

}
\item{out.fn}{NULL for no GDS output, or a file name}

\item{out.prec}{double or single precision for storage}

\item{out.compress}{the compression method for storing the GRM matrix in the GDS file}

\item{with.id}{if TRUE, the returned value with sample.id and sample.id}

\item{verbose}{if TRUE, show information}

\item{...}{passing to other SNP filtering parameters}

}
\details{The method estimates the z-score of each locus/allele relating to the multidimensional ancestry spaces.
If there are \strong{n} samples, there will be \strong{n X n} ancestry genetic trajectories, with an eigen decomposition, producting \strong{n} dimensional spaces that represent the common ancestry maps.
This method was conceived combining the idea of KLFDAPT (Qin, 2021) \url{https://xinghuq.github.io/KLFDAPC/articles/Genome_scan_KLFDAPC.html} and IBD-based genome scan (Albrechtsen et al., 2010). With an eigenvector decomposition of IBD (Zheng & Weir 2016), we can estimate the population ancestry propotion. It is competitive to pcadapt (Luu, 2016), as it considers \strong{n} latent genetic spaces (which is different from pcadapt that chooses \strong{k} components from \strong{p} eigenverctors).

}
\value{ A GenomeAdapt class, containing the loci z-scores and the eigen analysis results of IBD representing the genetic structure.
 \item{zscores}{The locus Z-scores relating to n latent ancestry genetic spaces, n is equal to the number of individuals}
  \item{eig}{Eigen analysis of IBD represent ancestry or population genetic structure}
  \item{chr }{Chromosomes }

}
\references{

Qin, X., Chiang, C. W., & Gaggiotti, O. E. (2021). Kernel Local Fisher Discriminant Analysis of Principal Components (KLFDAPC) significantly improves the accuracy of predicting geographic origin of individuals.bioRxiv.

Zheng, X., & Weir, B. S. (2016). Eigenanalysis of SNP data with an identity by descent interpretation. Theoretical population biology, 107, 65-76.

Albrechtsen, A., Moltke, I., & Nielsen, R. (2010). Natural selection and the distribution of identity-by-descent in the human genome. Genetics, 186(1), 295-308.

Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G. (2016). Detecting genomic signatures of natural selection with principal component analysis: application to the 1000 genomes data. Molecular biology and evolution, 33(4), 1082-1093.

}
\author{qin.xinghu@163.com

}
\examples{
### using an example dataset (HapMap)to conduct the genome scan

HapmapScan=GenomeAdapt.gds(genfile = SNPRelate::snpgdsExampleFileName(),method="EIGMIX",
num.thread = 1L, autosome.only=TRUE, remove.monosnp=TRUE, maf=0.01, missing.rate=0.1)

}

\keyword{GenomeAdapt}
