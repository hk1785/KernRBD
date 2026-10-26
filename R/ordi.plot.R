ordi.plot <-
function(K, Treatment, dim = 2, pcs = NULL, cols = NULL, pchs = NULL, main = NULL, tag = NULL,
         ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "topright", cex = 1, ...) {

  # Input checks

  K <- as.matrix(K)
  N <- nrow(K)
  if (!is.numeric(K) || ncol(K) != N || length(Treatment) != N) {
    stop("K should be a numeric N x N matrix, where N is the length of Treatment.")
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

  ind <- !is.na(Treatment)
  K <- K[ind, ind]
  Treatment <- factor(Treatment[ind])
  N <- nrow(K)
  grps <- levels(Treatment)
  m <- length(grps)

  if (is.null(cols)) {
    if (m <= 4) {
      cols <- c("green2", "blue2", "purple2", "red2")[1:m]
    } else {
      cols <- grDevices::hcl.colors(m, palette = "Dark 3")
    }
  }
  if (is.null(pchs)) {
    pchs <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 3, 4, 8)[(0:(m - 1)) %% 12 + 1]
  }
  cols <- rep(cols, length.out = m)
  pchs <- rep(pchs, length.out = m)

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
    stop("pcs exceed the number of principal coordinates.")
  }
  colnames(pc.mat) <- paste("PC", 1:ncol(pc.mat), sep = "")
  names(prop) <- colnames(pc.mat)
  labs <- paste("PC ", pcs, " (", sprintf("%.1f", prop[pcs] * 100), "%)", sep = "")

  col.obs <- cols[as.integer(Treatment)]
  pch.obs <- pchs[as.integer(Treatment)]

  if (dim == 2) {

    x <- pc.mat[, pcs[1]]
    y <- pc.mat[, pcs[2]]
    xr <- range(x)
    yr <- range(y)
    graphics::plot(x, y, type = "n", xlab = labs[1], ylab = labs[2], main = main,
                   xlim = xr + c(-1, 1) * diff(xr)/10, ylim = yr + c(-1, 1) * diff(yr)/10,
                   las = 1, mgp = c(2.5, 0.7, 0), cex.lab = cex, cex.axis = 0.9 * cex, cex.main = 1.1 * cex, ...)
    graphics::abline(h = 0, v = 0, col = "grey85", lty = 2)

    for (g in 1:m) {
      ind.g <- which(as.integer(Treatment) == g)
      cx <- mean(x[ind.g])
      cy <- mean(y[ind.g])
      if (star) {
        graphics::segments(cx, cy, x[ind.g], y[ind.g], col = grDevices::adjustcolor(cols[g], alpha.f = 0.4))
      }
      if (ellipse && length(ind.g) > 2) {
        S <- stats::cov.wt(cbind(x[ind.g], y[ind.g]), method = "ML")$cov
        e <- eigen(S, symmetric = TRUE)
        th <- seq(0, 2 * pi, length.out = 200)
        el <- t(c(cx, cy) + 1.5 * e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% rbind(cos(th), sin(th)))
        graphics::lines(el, col = cols[g], lwd = 1.5)
      }
    }
    graphics::points(x, y, col = col.obs, pch = pch.obs, cex = cex)
    if (legend) {
      graphics::legend(legend.pos, legend = grps, col = cols, pch = pchs, bty = "n", cex = 0.9 * cex)
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
    rgl::plot3d(x = pc.mat[, pcs[1]], y = pc.mat[, pcs[2]], z = pc.mat[, pcs[3]], col = col.obs, type = "s",
                xlab = labs[1], ylab = labs[2], zlab = labs[3], main = main, ...)
    if (legend) {
      rgl::legend3d(legend.pos, legend = grps, col = cols, pch = 15, bty = "n", cex = 1.3 * cex)
    }

  }

  invisible(list(pcs = pc.mat, prop = prop))
}
