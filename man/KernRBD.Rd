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
The number of resamples (Default: 5000).
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
  
  require(ade4)
  require(gtools)
  require(phyloseq)
  require(rgl)
  require(KernRBD)
  
  ###################
  # Perform KernRBD #
  ###################
  
  set.seed(521)
  
  s.time <- proc.time()
  out <- KernRBD(Ks = Ks, Block.IDs = sam.dat$Monkey.ID, Treatment = sam.dat$Treatment, n.res = 5000) 
  e.time <- proc.time()
  e.time - s.time
  
  # Global Testing
  
  round(out$global.out$itembyitem.pvals, 3)
  round(out$global.out$minP.pval, 3)
  
  # Pairwise Comparisons
  
  round(out$pairwise.out$itembyitem.adj.pvals, 3)
  round(out$pairwise.out$minP.adj.pvals, 3)
  
  #######################
  # Ordination Plotting #
  #######################
  
  ##### First kernel (Jaccard) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[1]]) # First kernel (Jaccard)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "Jaccard", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Second kernel (Bray-Curtis) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[2]]) # Second kernel (Bray-Curtis)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "Bray-Curtis", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Third kernel (UUniFrac) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[3]]) # Third kernel (UUniFrac)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "UUniFrac", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Fourth kernel (GUniFrac 0.25) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[4]]) # Fourth kernel (GUniFrac 0.25)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "GUniFrac (0.25)", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Fifth kernel (GUniFrac 0.5) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[5]]) # Fifth kernel (GUniFrac 0.5)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "GUniFrac (0.5)", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Sixth kernel (GUniFrac 0.75) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[6]]) # Sixth kernel (GUniFrac 0.75)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "GUniFrac (0.75)", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
  
  ### Seventh kernel (WUniFrac) 
  
  # Color labels for the treatments: "Baseline", "PRF 1", "PRF 2", "PRF 3"
  
  cols <- c("green2", "blue2", "purple2", "red2")
  
  # Singular value decomposition
  
  svd.out <- svd(out$Ks[[7]]) # Seventh kernel (WUniFrac)
  pcs <- (svd.out$u %*% sqrt(diag(svd.out$d)))[,svd.out$d > 0]
  eivals <- svd.out$d[svd.out$d > 0]
  
  # 2D plot
  
  plot(pcs[,1], pcs[,2], 
       xlab = paste("PC 1 (", round((eivals[1]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       ylab = paste("PC 2 (", round((eivals[2]/sum(eivals)) * 100, 1), "%)", sep = ""), 
       main = "WUniFrac", 
       col = cols, mgp = c(2.5, 1, 0), pch = 21, frame = TRUE,
       xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
       ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
       cex = 1, cex.main = 1.2, cex.lab = 1)
  
  s.class(as.data.frame(pcs[,c(1, 2)]), fac = as.factor(out$treatment), cpoint = 0, col = cols, add.plot = TRUE)
  
  # 3D plot
  
  plot3d( 
    x = pcs[,2], y = pcs[,3], z = pcs[,1], col = cols, 
    type = "s", 
    zlab = paste("PC 1 (", round(eivals[1]/sum(eivals) * 100, 1), "%)", sep = ""),
    xlab = paste("PC 2 (", round(eivals[2]/sum(eivals) * 100, 1), "%)", sep = ""),
    ylab = paste("PC 3 (", round(eivals[3]/sum(eivals) * 100, 1), "%)", sep = ""),
    cex.axis = 0.5)
  
  rglwidget()
}
