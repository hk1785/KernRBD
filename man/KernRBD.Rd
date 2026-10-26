\name{KernRBD}
\alias{KernRBD}
\title{
KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons
}
\description{
This function conducts KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons, for each input kernel and the omnibus test across multiple input kernels.
}
\usage{
KernRBD(Ks, Block.IDs, Treatment, n.res = 5000, n.cores = 1)
}
\arguments{
  \item{Ks}{
A list of kernel (e.g., similarity) matrices. A single kernel matrix is also allowed.
  }
  \item{Block.IDs}{
A vector of block (e.g., subject) IDs.
  }
  \item{Treatment}{
A vector of treatment labels.
  }
  \item{n.res}{
The number of resamples (i.e., within-block rearrangements) (Default: 5000). The smallest attainable P-value is 1/(n.res + 1); hence, a larger number (e.g., n.res = 10000 or more) is recommended for stringent significance levels after the adjustment for multiple comparisons.
  }
  \item{n.cores}{
The number of cores for parallel computation (Default: 1). This requires OpenMP; otherwise, a single core is used. The results do not depend on n.cores.
  }
}
\details{
The P-values for global testing are calculated by rearranging the treatment labels within each block, and the P-values for pairwise comparisons are calculated by rearranging only the two treatment labels of interest within each block. The adjusted P-values for pairwise comparisons are calculated using the step-down minP procedure. The omnibus test is based on the minimum P-value test statistic across the input kernels, where the same rearrangements are used for all input kernels. The kernels are centered, symmetrized and normalized by their traces, which does not change the P-values. For the details, see Koh (In Review).

The test statistic values are computed in C++ through the R package Rcpp, where the rearrangements are identical to those by sample() of base R for the same random seed, and they can be computed in parallel across multiple cores through OpenMP. The computational cost is O(|K|(1 + q)RN^2), where |K| is the number of input kernels, q is the number of pairwise comparisons, R is the number of resamples and N is the number of observations.
}
\value{
Ks: A list of kernel (similarity) matrices after quality controls and reorganizations.

block.IDs: A vector of block (e.g., subject) IDs after quality controls and reorganizations.

treatment: A vector of treatment labels after quality controls and reorganizations.

global.out: A list of P-values for global testing (itembyitem.pvals: P-values for each kernel; minP.pval: P-value for the omnibus test).

pairwise.out: A list of P-values and adjusted P-values for pairwise comparisons (itembyitem.pvals and itembyitem.adj.pvals: P-values and adjusted P-values for each kernel; minP.pvals and minP.adj.pvals: P-values and adjusted P-values for the omnibus test). This is provided for more than two treatments.
}
\references{
Koh H. A general kernel machine comparative analysis framework for randomized block designs. (In Review)
}
\author{
Hyunwook Koh
}
\examples{
data(Ds)
Ks <- lapply(Ds, function(d) D.to.K(d))

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(biom)
  sam.dat <- phyloseq::sample_data(biom)
  set.seed(521)
  out <- KernRBD(Ks = Ks, Block.IDs = sam.dat$Monkey.ID, Treatment = sam.dat$Treatment, n.res = 199)
  out$global.out
  out$pairwise.out$minP.adj.pvals
}
}
