\name{D.to.K}
\alias{D.to.K}
\title{
Distance to kernel
}
\description{
This function converts distance matrix to kernel matrix, i.e., K = -(1/2)(I - 11'/N)D^2(I - 11'/N), where D^2 is the matrix of squared distances.
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
Koh H. A general kernel machine comparative analysis framework for randomized block designs: Applications to human microbiome studies. (In Review)
}
\author{
Hyunwook Koh
}
\examples{
data(Ds)
Ks <- lapply(Ds, function(d) D.to.K(d))
}
