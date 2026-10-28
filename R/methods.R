# Whether a matrix looks like a distance matrix (zero diagonal and nonnegative entries) rather than a kernel matrix

is.distance.like <- function(M) {
  tol <- sqrt(.Machine$double.eps) * max(1, max(abs(M)))
  all(abs(diag(M)) <= tol) && all(M >= -tol) && any(M > tol)
}

# P-values formatted with a fixed number of digits, where small P-values are shown as "<0.001" (digits = 3)

fmt.pvals <- function(p, digits = 3) {
  out <- formatC(p, format = "f", digits = digits)
  small <- !is.na(p) & p < 10^(-digits)
  out[small] <- paste("<", formatC(10^(-digits), format = "f", digits = digits), sep = "")
  out[is.na(p)] <- "-"
  out
}

summary.KernRBD <- function(object, adjusted = TRUE, ...) {
  kern.names <- names(object$global.out$itembyitem.pvals)
  tab <- data.frame(Global = c(object$global.out$itembyitem.pvals, object$global.out$minP.pval),
                    row.names = c(kern.names, "Omnibus"), check.names = FALSE)
  if (!is.null(object$pairwise.out)) {
    pw <- object$pairwise.out
    P <- if (adjusted) rbind(pw$itembyitem.adj.pvals, Omnibus = pw$minP.adj.pvals) else rbind(pw$itembyitem.pvals, Omnibus = pw$minP.pvals)
    pairs <- attr(object, "pairs")
    if (is.null(pairs)) {
      pairs <- t(utils::combn(as.character(sort(unique(object$treatment))), 2))
    }
    colnames(P) <- paste(pairs[, 1], "vs.", pairs[, 2])
    tab <- cbind(tab, as.data.frame(P, check.names = FALSE))
  }
  attr(tab, "adjusted") <- adjusted
  class(tab) <- c("summary.KernRBD", "data.frame")
  tab
}

print.summary.KernRBD <- function(x, digits = 3, ...) {
  pairwise <- ncol(x) > 1
  cat(if (!pairwise) "P-values for global testing" else if (isTRUE(attr(x, "adjusted")))
        "P-values for global testing and adjusted P-values (step-down minP) for pairwise comparisons"
      else "P-values for global testing and unadjusted P-values for pairwise comparisons", "\n\n", sep = "")
  out <- as.matrix(as.data.frame(lapply(unclass(x), fmt.pvals, digits = digits), check.names = FALSE))
  dimnames(out) <- list(rownames(x), colnames(x))
  print(out, quote = FALSE, right = TRUE)
  invisible(x)
}

print.KernRBD <- function(x, digits = 3, ...) {
  treat.levels <- sort(unique(x$treatment))
  kern.names <- names(x$global.out$itembyitem.pvals)
  cat("KernRBD: Kernel machine comparative analysis for randomized block designs\n\n")
  cat("Observations: ", length(x$treatment), " (", length(unique(x$block.IDs)), " blocks)", sep = "")
  if (!is.null(attr(x, "n.removed")) && attr(x, "n.removed") > 0) {
    cat("; ", attr(x, "n.removed"), " removed for missing values", sep = "")
  }
  cat("\nTreatments: ", length(treat.levels), " (", paste(treat.levels, collapse = ", "), ")\n", sep = "")
  cat("Kernels: ", length(kern.names), " (", paste(kern.names, collapse = ", "), ")\n", sep = "")
  if (!is.null(attr(x, "n.res"))) {
    cat("Resamples: ", attr(x, "n.res"), "\n", sep = "")
  }
  if (!is.null(x$pairwise.out)) {
    cat("Pairwise comparisons: ", if (is.null(attr(x, "ref"))) "all pairs" else paste(attr(x, "ref"), "vs. each of the other treatments"), "\n", sep = "")
  }
  cat("\n")
  print(summary(x), digits = digits)
  invisible(x)
}
