# KernRBD

Title: Kernel Machine Comparative Analysis for Randomized Block Designs

Version: 1.2

Date: 2026-10-28

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package provides facilities for a general kernel machine comparative analysis framework for randomized block designs, named as KernRBD, to investigate the effects of treatments (e.g., medical treatment, environmental exposure) on the underlying variants. KernRBD is unique in its range of functionalities, including the computation of P-values for global testing and adjusted P-values for pairwise comparisons, as well as the visual representation using an ordination plot in a reduced coordinate space. KernRBD is also (i) practical requiring only a kernel as input without the need to know of its underlying real variants, and (ii) robustly valid based on a resampling scheme not requiring the assumption of normality to be satisfied. Furthermore, the omnibus test of KernRBD based on the minimum P-value test statistic enables a unified and powerful significance testing across multiple input kernels.

License: GPL-3

## Reference

* Koh H. A general kernel machine comparative analysis framework for randomized block designs: Applications to human microbiome studies. (_In Review_)

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/KernRBD/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Software Dependencies

KernRBD is implemented in R (>= 3.6.0) and C++, and requires only the base R packages (stats, graphics, grDevices and utils) and the R package Rcpp. The R package rgl is optional for three-dimensional ordination plots, and the R package phyloseq is optional for the example data. The installation from GitHub requires a C++ compiler (Rtools for Windows; Xcode Command Line Tools for macOS), and the multicore computation requires OpenMP.

## Prerequisites

Rcpp
```
install.packages("Rcpp")
```
rgl (optional)
```
install.packages("rgl")
```
phyloseq (optional)
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
This function conducts KernRBD to compute P-values for global testing and adjusted P-values for pairwise comparisons, for each input kernel and the omnibus test across multiple input kernels.

### Syntax
```
KernRBD(Ks, Block.IDs, Treatment, n.res = 5000, n.cores = 1)

print(out, digits = 3)
summary(out, adjusted = TRUE)
```

### Arguments
* _Ks_ - A list of kernel (e.g., similarity) matrices. A single kernel matrix is also allowed.
* _Block.IDs_ - A vector of block (e.g., subject) IDs.
* _Treatment_ - A vector of treatment labels.
* _n.res_ - The number of resamples (Default: 5000).
* _n.cores_ - The number of cores for parallel computation (Default: 1).

### Values
* _Ks_ - A list of kernel (similarity) matrices after quality controls and reorganizations.
* _block.IDs_ - A vector of block (e.g., subject) IDs after quality controls and reorganizations.
* _treatment_ - A vector of treatment labels after quality controls and reorganizations.
* _global.out_ - A list of P-values for global testing.
* _pairwise.out_ - A list of P-values and adjusted P-values for pairwise comparisons.

The output (_out_) is printed as a summary of the design and a table of the P-values for global testing and adjusted P-values for pairwise comparisons, for each kernel and the omnibus test (Omnibus). The same table is returned by summary(out) as a data frame (adjusted = FALSE for the unadjusted P-values for pairwise comparisons), which can be saved by write.csv().

## :mag: D.to.K

### Description
This function converts a distance matrix (or each distance matrix in a list) to a kernel matrix.

### Syntax
```
D.to.K(D)
```

### Arguments
* _D_ - A distance matrix, a dist object, or a list of them.

### Values
* A kernel matrix (or a list of kernel matrices).

## :mag: ordi.plot

### Description
This function draws a two- or three-dimensional ordination plot using the principal coordinates of a kernel.

### Syntax
```
ordi.plot(K, Treatment, dim = 2, pcs = NULL, cols = NULL, pchs = NULL, main = NULL, tag = NULL,
          ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "auto", legend.title = NULL, cex = 1, ...)
```

### Arguments
* _K_ - A kernel (e.g., similarity) matrix.
* _Treatment_ - A vector of treatment labels.
* _dim_ - The number of dimensions: 2 or 3 (Default: 2). The three-dimensional plot requires the R package rgl.
* _pcs_ - The principal coordinates to plot (Default: c(1, 2) for dim = 2 and c(1, 2, 3) for dim = 3).
* _cols_ - The colors for the treatments (Default: the colorblind-safe Okabe-Ito palette).
* _pchs_ - The symbols for the treatments (Default: filled symbols).
* _main_ - The title of the plot.
* _tag_ - The panel label (e.g., "A").
* _ellipse_ - Whether to draw an ellipse (1.5 standard deviations) for each treatment (Default: TRUE). A positive number sets the number of standard deviations (e.g., ellipse = 2).
* _star_ - Whether to connect the observations to their treatment centroid (Default: TRUE).
* _legend_ - Whether to draw a legend (Default: TRUE).
* _legend.pos_ - The position of the legend (Default: "auto", the corner overlapping the fewest observations), or a keyword such as "topleft".
* _legend.title_ - The title of the legend (e.g., "Treatment").
* _cex_ - The scaling factor for the symbols and texts (Default: 1).
* _..._ - Other graphical parameters passed to plot (dim = 2; e.g., xlim, ylim and asp) or rgl::plot3d (dim = 3).

### Values
* _pcs_ - A matrix of the principal coordinates.
* _prop_ - A vector of the proportions of total variance explained.

---------------------------------------------------------------------------------------------------------------------------------------

## Example

Import requisite R packages
```
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
Convert pairwise distance matrix to pairwise kernel matrix. The seven distances are Jaccard, Bray-Curtis, unweighted UniFrac, generalized UniFrac (0.25), generalized UniFrac (0.5), generalized UniFrac (0.75) and weighted UniFrac distances.
```
data(Ds)

Ks <- D.to.K(Ds)
names(Ks) <- c("Jaccard", "Bray-Curtis", "U.UniFrac", "G.UniFrac(0.25)", "G.UniFrac(0.5)", "G.UniFrac(0.75)", "W.UniFrac")
```
Perform KernRBD.
```
set.seed(521)
out <- KernRBD(Ks = Ks, Block.IDs = sam.dat$Monkey.ID, Treatment = sam.dat$Treatment, n.res = 5000)

# P-values for global testing and adjusted P-values for pairwise comparisons

out

# The same table as a data frame (e.g., to save as a CSV file)

tab <- summary(out)
write.csv(tab, "KernRBD_Results.csv")
```
Draw ordination plots.
```
# Two-dimensional plots (saved as a PDF file)

pdf("KernRBD_Ordination.pdf", width = 16, height = 8)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (i in 1:7) {
  ordi.plot(Ks[[i]], Treatment = sam.dat$Treatment, main = names(Ks)[i], tag = LETTERS[i], legend = (i == 1))
}
dev.off()

# Three-dimensional plot (First kernel: Jaccard)

ordi.plot(Ks[[1]], Treatment = sam.dat$Treatment, dim = 3)
```

---------------------------------------------------------------------------------------------------------------------------------------

## Recommended Settings

* _n.res_ - The smallest attainable P-value is 1/(n.res + 1). The default (n.res = 5000) is sufficient in most cases, but a larger number (e.g., n.res = 10000 or more) is recommended for stringent significance levels after the adjustment for multiple comparisons. A random seed (e.g., set.seed(521)) is recommended for reproducibility.
* _n.cores_ - Multiple cores (e.g., n.cores = parallel::detectCores() - 1) are recommended for large data. The results do not depend on n.cores for the same random seed.
* _Ks_ - The result can vary depending on the choice of kernel. Hence, we recommend multiple input kernels that represent distinct perspectives, together with the omnibus test. For example, for microbiome data, the Jaccard and unweighted UniFrac kernels (presence/absence), the Bray-Curtis and weighted UniFrac kernels (abundance), and the generalized UniFrac kernels (in-between), either incorporating or ignoring phylogenetic tree information. For other high-dimensional data (e.g., metabolomics, transcriptomics, proteomics), linear or Gaussian kernels of standardized features can be used. Adding many uninformative kernels can reduce power.

## Execution Time

The computational cost is O(|K|(1 + q)RN^2), where |K| is the number of input kernels, q is the number of pairwise comparisons, R is the number of resamples (n.res) and N is the number of observations. The test statistic values are computed in C++, and can be computed in parallel across multiple cores.

The runtimes of KernRBD per dataset (seconds) are summarized below for m treatments and N observations. The runtimes are the means (standard deviations) over 100 simulated datasets for the global testing and pairwise comparisons for seven input kernels and the omnibus test, with n.res = 5000 on a single core of a Mac Studio (Apple M1 Ultra; 64 GB of memory; macOS 26.3.1), excluding the computation of the input kernels.

| m | N | Runtime (seconds) |
|:---:|:---:|:---:|
| 3 | 75 | 0.13 (0.01) |
| 3 | 150 | 0.45 (0.01) |
| 4 | 100 | 0.24 (0.01) |
| 4 | 200 | 0.82 (0.03) |
| 5 | 125 | 0.35 (0.02) |
| 5 | 250 | 1.25 (0.06) |

## Version History

* Version 1.2 (2026-10-28) - Summary table of the P-values by print() and summary(), messages for removed observations and warnings for distance matrices given as kernels, D.to.K for a list of distance matrices, and the improved ordination plot (colorblind-safe colors, filled symbols, automatic legend position, legend title, and user-specified ellipse size and axis limits).
* Version 1.1 (2026-10-26) - Faster computation in C++ (Rcpp) with multicore computation (OpenMP), input checks, the ordination plot function (ordi.plot), and the P-values exactly following the formulas in the reference (e.g., ties counted in the reference set, and the step-down minP procedure based on the null P-values).
