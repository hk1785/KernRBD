# KernRBD

Title: Kernel Machine Comparative Analysis for Randomized Block Designs

Version: 1.0

Date: 2025-02-06

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package provides facilities for a general kernel machine comparative analysis framework for randomized block designs, named as KernRBD, to investigate the effects of treatments (e.g., medical treatment, environmental exposure) on the underlying variants. KernRBD is unique in its range of functionalities, including the computation of P-values for global testing and adjusted P-values for pairwise comparisons, as well as the visual representation using an ordination plot in a reduced coordinate space. KernRBD is also (i) practical requiring only a kernel as input without the need to know of its underlying real variants, and (ii) robustly valid based on a resampling scheme not requiring the assumption of normality to be satisfied. Furthermore, the omnibus test of KernRBD based on the minimum P-value test statistic enables a unified and powerful significance testing across multiple input kernels.

NeedsCompilation: no

Depends: R(>= 4.4.1), ade4, gtools, phyloseq, rgl

License: GPL v3.0

NeedsCompilation: no

## Reference

* Koh H. A general kernel machine comparative analysis framework for randomized block designs. (_In Review_)

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/KernRBD/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Prerequites

ade4, gtools, rgl
```
install.packages(c("ade4", "gtools", "rgl"))
```
phyloseq
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
```

## Installation

```
library(devtools)
install_github("hk1785/KernRBD", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------

## :mag: KernRBD

### Description
This function conducts KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons.

### Syntax
```
KernRBD(Ks, Block.IDs, Treatment, n.res = 5000)
```

### Inputs
* _Ks_ - A list of kernel (e.g., similarity) matrices.
* _Block.IDs_ - A vector of block (e.g., subject) IDs.
* _Treatment_ - A vector of treatment labels.
* _n.res_ - The number of resamples (Default: 5000).

### Outputs
* _Ks_ - A list of kernel (similarity) matrices after quality controls and reorganizations.
* _block.IDs_ - A vector of block (e.g., subject) IDs after quality controls and reorganizations.
* _treatment_ - A vector of treatment labels after quality controls and reorganizations.
* _global.out_ - A list of P-values for global testing.
* _pairwise.out_ - A list of P-values and adjusted P-values for pairwise comparisons.

### Example
Import requisite R packages
```
library(ade4)
library(gtools)
library(phyloseq)
library(rgl)
library(KernRBD)
```
Example Data: Gut microbiome data to assess the effect of periodically restricted feeding (PRF) on gut microbiome profiles (Yanai et al., 2024, Nat Commun)
```
data(biom)

otu.tab <- otu_table(biom)
tax.tab <- tax_table(biom)
sam.dat <- sample_data(biom)
tree <- phy_tree(biom)
```
Convert pairwise distance matrix to pairwise kernel matrix.
```
data(Ds)

Ks <- lapply(Ds, function(d) D.to.K(d))
```
Perform KernRBD.
```
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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
     col = cols, mgp = c(2.5, 1, 0), pch = 1:length(cols), frame = TRUE,
     xlim = c(min(pcs[,1]) - (max(pcs[,1]) - min(pcs[,1]))/10, max(pcs[,1]) + (max(pcs[,1]) - min(pcs[,1]))/10), 
     ylim = c(min(pcs[,2]) - (max(pcs[,2]) - min(pcs[,2]))/10, max(pcs[,2]) + (max(pcs[,2]) - min(pcs[,2]))/10), 
     cex = 1, cex.main = 1.2, cex.lab = 1)

grid(10, 10, lwd = 0.5)

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
```
