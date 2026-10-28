\name{KernRBD}
\alias{KernRBD}
\alias{print.KernRBD}
\alias{summary.KernRBD}
\alias{print.summary.KernRBD}
\title{
KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons
}
\description{
This function conducts KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons, for each input kernel and the omnibus test across multiple input kernels.
}
\usage{
KernRBD(Ks, Block.IDs, Treatment, ref = NULL, n.res = 5000, n.cores = 1)

\method{print}{KernRBD}(x, digits = 3, ...)

\method{summary}{KernRBD}(object, adjusted = TRUE, ...)
}
\arguments{
  \item{Ks}{
A list of kernel (e.g., similarity) matrices. A single kernel matrix is also allowed.
  }
  \item{Block.IDs}{
A vector of block (e.g., subject) IDs.
  }
  \item{Treatment}{
A vector of treatment labels (e.g., numbers, characters or factor levels), each of which is treated as a distinct treatment.
  }
  \item{ref}{
The reference treatment label (e.g., baseline) for pairwise comparisons (Default: NULL, i.e., all possible pairwise comparisons). If specified, only the pairwise comparisons between the reference treatment and each of the other treatments are conducted, and the adjusted P-values are calculated for these comparisons.
  }
  \item{n.res}{
The number of resamples (i.e., within-block rearrangements) (Default: 5000). The smallest attainable P-value is 1/(n.res + 1); hence, a larger number (e.g., n.res = 10000 or more) is recommended for stringent significance levels after the adjustment for multiple comparisons.
  }
  \item{n.cores}{
The number of cores for parallel computation (Default: 1). This requires OpenMP; otherwise, a single core is used. The results do not depend on n.cores.
  }
  \item{x, object}{
An output of KernRBD.
  }
  \item{digits}{
The number of decimal places for the P-values (Default: 3), where smaller P-values are shown as, e.g., <0.001.
  }
  \item{adjusted}{
Whether to report the adjusted P-values (TRUE) or the unadjusted P-values (FALSE) for pairwise comparisons (Default: TRUE).
  }
  \item{...}{
Not used.
  }
}
\details{
The P-values for global testing are calculated by rearranging the treatment labels within each block, and the P-values for pairwise comparisons are calculated by rearranging only the two treatment labels of interest within each block. The adjusted P-values for pairwise comparisons are calculated using the step-down minP procedure. The omnibus test is based on the minimum P-value test statistic across the input kernels, where the same rearrangements are used for all input kernels. The kernels are centered, symmetrized and normalized by their traces, which does not change the P-values. For the details, see Koh (In Review).

The output is printed as a summary of the design and a table of the P-values for global testing and adjusted P-values for pairwise comparisons, for each kernel (rows) and the omnibus test (last row). The same table is returned by summary() as a data frame, which can be saved, e.g., by write.csv(). Observations with missing Block.IDs or Treatment are removed with a message, and a warning is given if an input matrix looks like a distance matrix rather than a kernel matrix.

The test statistic values are computed in C++ through the R package Rcpp, where the rearrangements are identical to those by sample() of base R for the same random seed, and they can be computed in parallel across multiple cores through OpenMP. The computational cost is O(|K|(1 + q)RN^2), where |K| is the number of input kernels, q is the number of pairwise comparisons, R is the number of resamples and N is the number of observations.
}
\value{
Ks: A list of kernel (similarity) matrices after quality controls and reorganizations.

block.IDs: A vector of block (e.g., subject) IDs after quality controls and reorganizations.

treatment: A vector of treatment labels after quality controls and reorganizations.

global.out: A list of P-values for global testing (itembyitem.pvals: P-values for each kernel; minP.pval: P-value for the omnibus test).

pairwise.out: A list of P-values and adjusted P-values for pairwise comparisons (itembyitem.pvals and itembyitem.adj.pvals: P-values and adjusted P-values for each kernel; minP.pvals and minP.adj.pvals: P-values and adjusted P-values for the omnibus test). This is provided for more than two treatments.

summary() returns a data frame of the P-values for global testing (Global) and the (adjusted) P-values for pairwise comparisons (e.g., "A vs. B"), for each kernel and the omnibus test (Omnibus).
}
\references{
Koh H. A general kernel machine comparative analysis framework for randomized block designs: Applications to human microbiome studies. (In Review)
}
\author{
Hyunwook Koh
}
\examples{
data(Ds)
Ks <- D.to.K(Ds)

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(biom)
  sam.dat <- phyloseq::sample_data(biom)
  set.seed(521)
  out <- KernRBD(Ks = Ks, Block.IDs = sam.dat$Monkey.ID, Treatment = sam.dat$Treatment, n.res = 99)
  out.ref <- KernRBD(Ks = Ks, Block.IDs = sam.dat$Monkey.ID, Treatment = sam.dat$Treatment,
                     ref = 1, n.res = 99)
  out
  summary(out)
}
}
