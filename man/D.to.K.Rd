\name{D.to.K}
\alias{D.to.K}
\title{
Distance to kernel
}
\description{
This function converts distance matrix to kernel matrix.
}
\usage{
D.to.K(D)
}
\arguments{
  \item{D}{
Distance matrix
}
}
\value{
Kernel matrix
}
\references{
Koh H. A general kernel machine comparative analysis framework for randomized block designs. (In Review)
}
\author{
Hyunwook Koh
}
\examples{
Ks <- lapply(Ds, function(d) D.to.K(d))
}
