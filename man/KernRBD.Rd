\name{KernRBD}
\alias{KernRBD}
\title{
KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons
}
\description{
KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons
}
\usage{
KernRBD(Ks, Block.IDs, Treatment, n.res = 5000)
}
\arguments{
  \item{Ks}{
A list of kernel (e.g., similarity) matrices.
  }
  \item{Block.IDs}{
A vector of block (e.g., subject) IDs.
  }
  \item{Treatment}{
A vector of treatment labels.
  }
  \item{n.res}{
The number of resamples (Default: 5000). Using a larger number (e.g., n.res = 30000) is recommended.
  }
}

\value{
Ks: A list of kernel (similarity) matrices after quality controls and reorganizations.
block.IDs: A vector of block (e.g., subject) IDs after quality controls and reorganizations.
treatment: A vector of treatment labels after quality controls and reorganizations.
global.out: A list of P-values for global testing.
pairwise.out: A list of P-values and adjusted P-values for pairwise comparisons.
}

\references{
Koh H. A general kernel machine comparative analysis framework for randomized block designs. (In Review)
}

\author{
Hyunwook Koh
}

\examples{
}
