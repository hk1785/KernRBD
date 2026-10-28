ordi.plot <-
function(K, Treatment, dim = 2, pcs = NULL, cols = NULL, pchs = NULL, main = NULL, tag = NULL,
         ellipse = TRUE, star = TRUE, legend = TRUE, legend.pos = "auto", legend.title = NULL,
         cex = 1, cex.pt = 1.2, cex.text = 1.1, theta = 40, phi = 20, interactive = FALSE, ...) {

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
  s.pt <- cex * cex.pt          # Size of the symbols
  s.tx <- cex * cex.text        # Size of the texts

  ind <- !is.na(Treatment)
  obs.names <- rownames(K)[ind]
  K <- K[ind, ind, drop = FALSE]
  Treatment <- Treatment[ind]
  Treatment <- if (is.factor(Treatment)) droplevels(Treatment) else factor(Treatment)
  N <- nrow(K)
  grps <- levels(Treatment)
  m <- length(grps)

  # Colors (distinct also in grayscale and for color vision deficiencies, for up to seven treatments) and symbols
  # (filled for up to five treatments)

  if (is.null(cols)) {
    pal <- c("#0072B2", "#E69F00", "#009E73", "#882255", "#56B4E9", "#D55E00", "#000000")
    cols <- if (m <= 7) pal[1:m] else grDevices::hcl.colors(m, palette = "Dark 3")
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
                      mgp = c(2.5, 0.7, 0) * max(1, s.tx), cex.lab = s.tx, cex.axis = 0.9 * s.tx, cex.main = 1.1 * s.tx)
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
    graphics::points(x, y, pch = pchs[g.obs], col = borders[g.obs], bg = cols[g.obs], cex = s.pt, lwd = 0.8)
    leg.xy <- cbind(x, y)
    leg.ell <- do.call(rbind, ells)

  } else if (!interactive) {

    # Three-dimensional plot (perspective), where the observations are dropped to the floor of the box for depth

    xyz <- pc.mat[, pcs]
    lims <- apply(xyz, 2, function(v) range(v) + c(-1, 1) * diff(range(v))/20)
    persp.args <- list(x = lims[, 1], y = lims[, 2], z = matrix(lims[1, 3], 2, 2), zlim = lims[, 3], theta = theta, phi = phi,
                       r = 4, expand = 0.9, border = NA, col = NA, box = TRUE, ticktype = "detailed", nticks = 4,
                       xlab = "", ylab = "", zlab = "", cex.axis = 0.8 * s.tx, lwd = 0.8)
    dots <- list(...)
    persp.args[names(dots)] <- dots
    pmat <- do.call(graphics::persp, persp.args)
    xy <- function(v) grDevices::trans3d(v[1], v[2], v[3], pmat)

    # Floor grid, drop lines and observations (from the back to the front)

    gx <- pretty(lims[, 1], 4)
    gy <- pretty(lims[, 2], 4)
    for (v in gx[gx >= lims[1, 1] & gx <= lims[2, 1]]) graphics::lines(grDevices::trans3d(c(v, v), lims[, 2], lims[1, 3], pmat), col = "grey88")
    for (v in gy[gy >= lims[1, 2] & gy <= lims[2, 2]]) graphics::lines(grDevices::trans3d(lims[, 1], c(v, v), lims[1, 3], pmat), col = "grey88")
    top <- grDevices::trans3d(xyz[, 1], xyz[, 2], xyz[, 3], pmat)
    bot <- grDevices::trans3d(xyz[, 1], xyz[, 2], lims[1, 3], pmat)
    graphics::segments(bot$x, bot$y, top$x, top$y, col = grDevices::adjustcolor(cols[g.obs], alpha.f = 0.45))
    graphics::points(bot$x, bot$y, pch = 16, cex = 0.25 * s.pt, col = grDevices::adjustcolor(cols[g.obs], alpha.f = 0.45))
    o <- order((cbind(xyz, 1) %*% pmat)[, 4], decreasing = TRUE)
    graphics::points(top$x[o], top$y[o], pch = pchs[g.obs][o], bg = cols[g.obs][o], col = borders[g.obs][o], cex = s.pt, lwd = 0.8)

    # Axis labels outside the tick labels, parallel to the axes

    r <- lims[2, ] - lims[1, ]
    angle <- function(a, b) {
      d <- c(diff(graphics::grconvertX(c(xy(a)$x, xy(b)$x), "user", "inches")), diff(graphics::grconvertY(c(xy(a)$y, xy(b)$y), "user", "inches")))
      ang <- atan2(d[2], d[1]) * 180/pi
      if (ang > 90) ang - 180 else if (ang < -90) ang + 180 else ang
    }
    corners <- expand.grid(x = lims[, 1], y = lims[, 2])
    scr <- grDevices::trans3d(corners$x, corners$y, lims[1, 3], pmat)
    fr <- unlist(corners[which.min(scr$y), ])
    lf <- unlist(corners[which.min(scr$x), ])
    out <- function(v, k) v + k * ifelse(v == lims[1, 1:2], -1, 1) * r[1:2]
    graphics::text(xy(c(mean(lims[, 1]), out(fr, 0.32)[2], lims[1, 3])), labs[1], cex = s.tx, xpd = NA,
                   srt = angle(c(lims[1, 1], fr[2], lims[1, 3]), c(lims[2, 1], fr[2], lims[1, 3])))
    graphics::text(xy(c(out(fr, 0.32)[1], mean(lims[, 2]), lims[1, 3])), labs[2], cex = s.tx, xpd = NA,
                   srt = angle(c(fr[1], lims[1, 2], lims[1, 3]), c(fr[1], lims[2, 2], lims[1, 3])))
    graphics::text(xy(c(out(lf, 0.13), mean(lims[, 3]))), labs[3], cex = s.tx, xpd = NA, srt = 90)
    if (!is.null(main)) {
      graphics::title(main = main, line = 1, cex.main = 1.1 * s.tx)
    }
    leg.xy <- cbind(top$x, top$y)
    box.pts <- do.call(rbind, lapply(1:3, function(k) {     # points on the edges of the box (for the legend position)
      do.call(rbind, lapply(1:4, function(e) {
        ends <- as.matrix(expand.grid(1:2, 1:2))[e, ]
        v <- matrix(NA, 20, 3)
        v[, k] <- seq(lims[1, k], lims[2, k], length.out = 20)
        v[, -k] <- matrix(c(lims[ends[1], (1:3)[-k][1]], lims[ends[2], (1:3)[-k][2]]), 20, 2, byrow = TRUE)
        p <- grDevices::trans3d(v[, 1], v[, 2], v[, 3], pmat)
        cbind(p$x, p$y)
      }))
    }))
    leg.ell <- box.pts

  } else {

    # Interactive three-dimensional plot (rgl)

    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("the R package rgl is needed for interactive three-dimensional plots.")
    }
    if (rgl::cur3d() == 0) {
      rgl::open3d(windowRect = c(50, 50, 650, 650))
    }
    rgl::par3d(cex = s.tx)
    rgl::plot3d(x = pc.mat[, pcs[1]], y = pc.mat[, pcs[2]], z = pc.mat[, pcs[3]], col = cols[g.obs], type = "s", size = 3 * s.pt,
                xlab = labs[1], ylab = labs[2], zlab = labs[3], main = main, ...)
    if (legend) {
      rgl::legend3d(if (identical(legend.pos, "auto")) "topright" else legend.pos, legend = grps, col = cols, pch = 16,
                    bty = "n", cex = 1.2 * s.tx, title = legend.title)
    }

  }

  # Legend and panel label (two-dimensional and three-dimensional plots)

  if (dim == 2 || !interactive) {
    if (legend) {
      leg.args <- list(legend = grps, pch = pchs, col = borders, pt.bg = cols, pt.lwd = 0.8, pt.cex = s.pt, bty = "n",
                       cex = 0.9 * s.tx, title = legend.title, inset = 0.01)
      if (identical(legend.pos, "auto")) {
        legend.pos <- if (dim == 2) auto.legend.pos(leg.args, leg.xy, leg.ell) else auto.legend.pos(leg.args, leg.xy, leg.ell, top.only = TRUE)
      }
      do.call(graphics::legend, c(list(legend.pos), leg.args))
    }
    if (!is.null(tag)) {
      graphics::mtext(tag, side = 3, line = 1.5, adj = 0, at = graphics::grconvertX(0.02, "nfc", "user"), font = 2, cex = 1.3 * s.tx)
    }
  }

  invisible(list(pcs = pc.mat, prop = prop))
}

# Darker version of colors (for the borders of the filled symbols)

darken <- function(cols, f = 0.65) {
  rgb.mat <- grDevices::col2rgb(cols) * f
  grDevices::rgb(rgb.mat[1, ], rgb.mat[2, ], rgb.mat[3, ], maxColorValue = 255)
}

# Legend position among the four (or two top) corners that overlaps the fewest observations (and ellipses or box edges)

auto.legend.pos <- function(leg.args, xy, ell.xy = NULL, top.only = FALSE) {
  corners <- if (top.only) c("topright", "topleft") else c("topright", "topleft", "bottomright", "bottomleft")
  n.in <- function(pts, r) {
    if (is.null(pts)) return(0)
    sum(pts[, 1] >= r$left & pts[, 1] <= r$left + r$w & pts[, 2] <= r$top & pts[, 2] >= r$top - r$h)
  }
  score <- sapply(corners, function(p) {
    r <- do.call(graphics::legend, c(list(p), leg.args, list(plot = FALSE)))$rect
    n.in(xy, r) + n.in(ell.xy, r)/(if (top.only) 4 else 20)
  })
  corners[which.min(score)]
}
