\name{D.to.K}
\alias{D.to.K}
\title{
Distance to kernel
}
\description{
This function converts a distance matrix (or each distance matrix in a list) to a kernel matrix, i.e., K = -(1/2)(I - 11'/N)D^2(I - 11'/N), where D^2 is the matrix of squared distances.
}
\usage{
D.to.K(D)
}
\arguments{
  \item{D}{
A distance matrix, a dist object, or a list of them.
}
}
\value{
A kernel matrix with the same row and column names as D, or a list of kernel matrices for a list of distance matrices.
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
}
