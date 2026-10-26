\name{ordi.plot}
\alias{ordi.plot}
\title{
Ordination plot
}
\description{
This function draws a two- or three-dimensional ordination plot using the principal coordinates of a kernel, where the observations are distinguished using different colors (and symbols) for different treatments.
}
\usage{
ordi.plot(K, Treatment, dim = 2, pcs = NULL, cols = NULL, pchs = NULL, main = NULL,
          tag = NULL, ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "topright",
          cex = 1, ...)
}
\arguments{
  \item{K}{
A kernel (e.g., similarity) matrix.
  }
  \item{Treatment}{
A vector of treatment labels.
  }
  \item{dim}{
The number of dimensions: 2 or 3 (Default: 2). The three-dimensional plot requires the R package rgl.
  }
  \item{pcs}{
The principal coordinates to plot (Default: c(1, 2) for dim = 2 and c(1, 2, 3) for dim = 3).
  }
  \item{cols}{
The colors for the treatments (Default: "green2", "blue2", "purple2" and "red2" for up to four treatments).
  }
  \item{pchs}{
The symbols for the treatments in the two-dimensional plot.
  }
  \item{main}{
The title of the plot.
  }
  \item{tag}{
The panel label (e.g., "A") at the top-left corner of the two-dimensional plot.
  }
  \item{ellipse}{
Whether to draw an ellipse (1.5 standard deviations) for each treatment in the two-dimensional plot (Default: TRUE).
  }
  \item{star}{
Whether to connect the observations to their treatment centroid in the two-dimensional plot (Default: TRUE).
  }
  \item{legend}{
Whether to draw a legend (Default: TRUE).
  }
  \item{legend.pos}{
The position of the legend (Default: "topright").
  }
  \item{cex}{
The scaling factor for the symbols and texts (Default: 1).
  }
  \item{...}{
Other graphical parameters passed to plot (dim = 2) or rgl::plot3d (dim = 3).
  }
}
\details{
The principal coordinates are obtained through the singular value decomposition of the kernel matrix, centered and normalized by its trace as in KernRBD, and the proportions of total variance explained are shown in the axis labels.
}
\value{
pcs: A matrix of the principal coordinates.

prop: A vector of the proportions of total variance explained.
}
\references{
Koh H. A general kernel machine comparative analysis framework for randomized block designs: Applications to human microbiome studies. (In Review)
}
\author{
Hyunwook Koh
}
\examples{
data(Ds)
K <- D.to.K(Ds[[1]])

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(biom)
  sam.dat <- phyloseq::sample_data(biom)
  ordi.plot(K, Treatment = sam.dat$Treatment, main = "Jaccard")
}
}
