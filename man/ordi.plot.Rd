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
          tag = NULL, ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "auto",
          legend.title = NULL, cex = 1, cex.pt = 1.2, cex.text = 1.1, theta = 40,
          phi = 20, interactive = FALSE, ...)
}
\arguments{
  \item{K}{
A kernel (e.g., similarity) matrix.
  }
  \item{Treatment}{
A vector of treatment labels.
  }
  \item{dim}{
The number of dimensions: 2 or 3 (Default: 2).
  }
  \item{pcs}{
The principal coordinates to plot (Default: c(1, 2) for dim = 2 and c(1, 2, 3) for dim = 3).
  }
  \item{cols}{
The colors for the treatments (Default: blue, orange, green and wine for up to four treatments).
  }
  \item{pchs}{
The symbols for the treatments (Default: circle, triangle, square and diamond for up to four treatments).
  }
  \item{main}{
The title of the plot.
  }
  \item{tag}{
The panel label (e.g., "A") at the top-left corner of the plot.
  }
  \item{ellipse}{
Whether to draw an ellipse (1.5 standard deviations) for each treatment in the two-dimensional plot (Default: TRUE). A positive number sets the number of standard deviations (e.g., ellipse = 2).
  }
  \item{star}{
Whether to connect the observations to their treatment centroid in the two-dimensional plot (Default: TRUE).
  }
  \item{legend}{
Whether to draw a legend (Default: TRUE).
  }
  \item{legend.pos}{
The position of the legend (Default: "auto", the corner overlapping the fewest observations), or a keyword such as "topleft" or "bottomright".
  }
  \item{legend.title}{
The title of the legend (e.g., "Treatment").
  }
  \item{cex}{
The overall scaling factor for the symbols and texts (Default: 1).
  }
  \item{cex.pt}{
The size of the symbols (Default: 1.2).
  }
  \item{cex.text}{
The size of the texts, i.e., the axis labels, tick labels, title, legend and panel label (Default: 1.1).
  }
  \item{theta, phi}{
The viewing angles (degrees) of the three-dimensional plot (Default: 40 and 20).
  }
  \item{interactive}{
Whether to draw the three-dimensional plot to be rotated and zoomed using the R package rgl (Default: FALSE).
  }
  \item{...}{
Other graphical parameters passed to plot (dim = 2; e.g., xlim, ylim and asp), persp (dim = 3) or rgl::plot3d (dim = 3 and interactive = TRUE).
  }
}
\details{
The principal coordinates are obtained through the singular value decomposition of the kernel matrix, centered and normalized by its trace as in KernRBD, and the proportions of total variance explained are shown in the axis labels. The three-dimensional plot shows each observation with a vertical line to the floor of the box for depth. For a publication-quality figure, the plot can be drawn on a graphics device such as pdf() or png(), and multiple panels can be arranged by par(mfrow = ...) with the tag argument.
}
\value{
pcs: A matrix of the principal coordinates (rows: observations; columns: PC1, PC2, ...).

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
  ordi.plot(K, Treatment = sam.dat$Treatment, main = "Jaccard", legend.title = "Treatment")
  ordi.plot(K, Treatment = sam.dat$Treatment, dim = 3, main = "Jaccard", legend.title = "Treatment")
}
}
