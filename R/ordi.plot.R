ordi.plot <-
function(K, Treatment, dim = 2, pcs = NULL, cols = NULL, pchs = NULL, main = NULL, tag = NULL,
         ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "auto", legend.title = NULL, cex = 1, ...) {

  # Input checks

  K <- as.matrix(K)
  N <- nrow(K)
  if (!is.numeric(K) || ncol(K) != N || length(Treatment) != N) {
    stop("K should be a numeric N x N matrix, where N is the length of Treatment.")
  }
  if (any(is.na(K))) {
    stop("K should not contain missing values.")
  }
  if (is.distance.like(K)) {
    warning("K looks like a distance matrix (zero diagonal and nonnegative entries); ",
            "convert a distance matrix to a kernel matrix using D.to.K().", call. = FALSE)
  }
  if (!is.element(dim, c(2, 3))) {
    stop("dim should be 2 or 3.")
  }
  if (is.null(pcs)) {
    pcs <- 1:dim
  }
  if (length(pcs) != dim) {
    stop("the length of pcs should be equal to dim.")
  }
  if (length(ellipse) != 1 || is.na(ellipse) || !(is.logical(ellipse) || (is.numeric(ellipse) && ellipse > 0))) {
    stop("ellipse should be TRUE, FALSE or a positive number (the number of standard deviations).")
  }
  ell.sd <- if (isTRUE(ellipse)) 1.5 else if (is.numeric(ellipse)) ellipse else 0

  ind <- !is.na(Treatment)
  obs.names <- rownames(K)[ind]
  K <- K[ind, ind, drop = FALSE]
  Treatment <- Treatment[ind]
  Treatment <- if (is.factor(Treatment)) droplevels(Treatment) else factor(Treatment)
  N <- nrow(K)
  grps <- levels(Treatment)
  m <- length(grps)

  # Colors (Okabe-Ito colorblind-safe palette for up to seven treatments) and symbols (filled for up to five treatments)

  if (is.null(cols)) {
    okabe.ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")
    cols <- if (m <= 7) okabe.ito[1:m] else grDevices::hcl.colors(m, palette = "Dark 3")
  }
  if (is.null(pchs)) {
    pchs <- c(21, 24, 22, 23, 25, 1, 2, 0, 5, 6, 3, 4, 8)[(0:(m - 1)) %% 13 + 1]
  }
  cols <- rep(cols, length.out = m)
  pchs <- rep(pchs, length.out = m)
  filled <- is.element(pchs, 21:25)
  borders <- ifelse(filled, darken(cols), cols)

  # Principal coordinates through the singular value decomposition of the centered kernel
  # normalized by its trace (Eq. 11)

  C <- diag(N) - matrix(1/N, N, N)
  K <- C %*% K %*% C
  K <- K/sum(diag(K))
  svd.out <- svd(K)
  pos <- svd.out$d > sqrt(.Machine$double.eps) * max(svd.out$d)
  eivals <- svd.out$d[pos]
  pc.mat <- (svd.out$u[, pos, drop = FALSE] %*% diag(sqrt(eivals), nrow = length(eivals)))
  prop <- eivals/sum(eivals)
  if (max(pcs) > ncol(pc.mat)) {
    stop("pcs exceed the number of principal coordinates (", ncol(pc.mat), ").")
  }
  colnames(pc.mat) <- paste("PC", 1:ncol(pc.mat), sep = "")
  rownames(pc.mat) <- obs.names
  names(prop) <- colnames(pc.mat)
  labs <- paste("PC ", pcs, " (", sprintf("%.1f", prop[pcs] * 100), "%)", sep = "")

  g.obs <- as.integer(Treatment)

  if (dim == 2) {

    x <- pc.mat[, pcs[1]]
    y <- pc.mat[, pcs[2]]

    # Centroids and ellipses (ell.sd standard deviations) for each treatment

    cents <- t(sapply(1:m, function(g) c(mean(x[g.obs == g]), mean(y[g.obs == g]))))
    ells <- lapply(1:m, function(g) {
      ind.g <- which(g.obs == g)
      if (ell.sd == 0 || length(ind.g) < 3) return(NULL)
      S <- stats::cov.wt(cbind(x[ind.g], y[ind.g]), method = "ML")$cov
      e <- eigen(S, symmetric = TRUE)
      th <- seq(0, 2 * pi, length.out = 200)
      t(cents[g, ] + ell.sd * e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% rbind(cos(th), sin(th)))
    })

    # Plotting region covering the observations and ellipses; xlim and ylim can be given in ...

    dots <- list(...)
    all.xy <- rbind(cbind(x, y), do.call(rbind, ells))
    xr <- range(all.xy[, 1])
    yr <- range(all.xy[, 2])
    if (is.null(dots$xlim)) dots$xlim <- xr + c(-1, 1) * diff(xr)/20
    if (is.null(dots$ylim)) dots$ylim <- yr + c(-1, 1) * diff(yr)/20
    plot.args <- list(x = x, y = y, type = "n", xlab = labs[1], ylab = labs[2], main = main, las = 1,
                      mgp = c(2.5, 0.7, 0), cex.lab = cex, cex.axis = 0.9 * cex, cex.main = 1.1 * cex)
    plot.args[names(dots)] <- dots
    do.call(graphics::plot, plot.args)
    graphics::abline(h = 0, v = 0, col = "grey85", lty = 2)

    for (g in 1:m) {
      ind.g <- which(g.obs == g)
      if (!is.null(ells[[g]])) {
        graphics::polygon(ells[[g]], col = grDevices::adjustcolor(cols[g], alpha.f = 0.08), border = NA)
      }
      if (star) {
        graphics::segments(cents[g, 1], cents[g, 2], x[ind.g], y[ind.g], col = grDevices::adjustcolor(cols[g], alpha.f = 0.4))
      }
      if (!is.null(ells[[g]])) {
        graphics::lines(ells[[g]], col = cols[g], lwd = 1.5)
      }
    }
    graphics::points(x, y, pch = pchs[g.obs], col = borders[g.obs], bg = cols[g.obs], cex = cex, lwd = 0.8)

    if (legend) {
      leg.args <- list(legend = grps, pch = pchs, col = borders, pt.bg = cols, pt.lwd = 0.8, bty = "n",
                       cex = 0.9 * cex, title = legend.title, inset = 0.01)
      if (identical(legend.pos, "auto")) {
        legend.pos <- auto.legend.pos(leg.args, cbind(x, y), do.call(rbind, ells))
      }
      do.call(graphics::legend, c(list(legend.pos), leg.args))
    }
    if (!is.null(tag)) {
      graphics::mtext(tag, side = 3, line = 1.5, adj = 0, at = graphics::grconvertX(0.02, "nfc", "user"), font = 2, cex = 1.3 * cex)
    }

  } else {

    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("the R package rgl is needed for three-dimensional plots.")
    }
    if (rgl::cur3d() == 0) {
      rgl::open3d(windowRect = c(50, 50, 650, 650))
    }
    rgl::plot3d(x = pc.mat[, pcs[1]], y = pc.mat[, pcs[2]], z = pc.mat[, pcs[3]], col = cols[g.obs], type = "s",
                xlab = labs[1], ylab = labs[2], zlab = labs[3], main = main, ...)
    if (legend) {
      if (identical(legend.pos, "auto")) {
        legend.pos <- "topright"
      }
      rgl::legend3d(legend.pos, legend = grps, col = cols, pch = 16, bty = "n", cex = 1.3 * cex, title = legend.title)
    }

  }

  invisible(list(pcs = pc.mat, prop = prop))
}

# Darker version of colors (for the borders of the filled symbols)

darken <- function(cols, f = 0.65) {
  rgb.mat <- grDevices::col2rgb(cols) * f
  grDevices::rgb(rgb.mat[1, ], rgb.mat[2, ], rgb.mat[3, ], maxColorValue = 255)
}

# Legend position among the four corners that overlaps the fewest observations (and ellipses)

auto.legend.pos <- function(leg.args, xy, ell.xy = NULL) {
  corners <- c("topright", "topleft", "bottomright", "bottomleft")
  n.in <- function(pts, r) {
    if (is.null(pts)) return(0)
    sum(pts[, 1] >= r$left & pts[, 1] <= r$left + r$w & pts[, 2] <= r$top & pts[, 2] >= r$top - r$h)
  }
  score <- sapply(corners, function(p) {
    r <- do.call(graphics::legend, c(list(p), leg.args, list(plot = FALSE)))$rect
    n.in(xy, r) + n.in(ell.xy, r)/20
  })
  corners[which.min(score)]
}
