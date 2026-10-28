KernRBD <-
function(Ks, Block.IDs, Treatment, n.res = 5000, n.cores = 1) {

  # Input checks

  if (is.matrix(Ks)) {
    Ks <- list(Ks)
  }
  if (!is.list(Ks) || length(Ks) == 0) {
    stop("Ks should be a list of kernel matrices.")
  }
  if (length(Block.IDs) != length(Treatment)) {
    stop("Block.IDs and Treatment should have the same length.")
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- as.matrix(Ks[[i]])
    if (!is.numeric(Ks[[i]]) || nrow(Ks[[i]]) != ncol(Ks[[i]]) || nrow(Ks[[i]]) != length(Treatment)) {
      stop("each kernel matrix should be a numeric N x N matrix, where N is the length of Treatment.")
    }
    if (any(is.na(Ks[[i]]))) {
      stop("kernel matrices should not contain missing values.")
    }
    if (max(abs(Ks[[i]] - t(Ks[[i]]))) > sqrt(.Machine$double.eps) * max(1, max(abs(Ks[[i]])))) {
      stop("kernel matrices should be symmetric.")
    }
    if (is.distance.like(Ks[[i]])) {
      warning("Ks[[", i, "]] looks like a distance matrix (zero diagonal and nonnegative entries); ",
              "convert distance matrices to kernel matrices using D.to.K().", call. = FALSE)
    }
  }
  if (length(n.res) != 1 || !is.numeric(n.res) || is.na(n.res) || n.res < 1) {
    stop("n.res should be a positive integer.")
  }
  n.res <- as.integer(n.res)
  if (length(n.cores) != 1 || !is.numeric(n.cores) || is.na(n.cores) || n.cores < 1) {
    stop("n.cores should be a positive integer.")
  }
  n.cores <- as.integer(n.cores)
  if (n.cores > 1 && !kernrbd_openmp()) {
    warning("OpenMP is not available; a single core is used.")
    n.cores <- 1L
  }

  kern.names <- names(Ks)
  if (is.null(kern.names) || any(kern.names == "")) {
    kern.names <- paste("Kern.", 1:length(Ks), sep = "")
  }

  # Quality controls and reorganizations

  ind.na <- which(is.na(Block.IDs) | is.na(Treatment))
  if (length(ind.na) > 0) {
    message(length(ind.na), " observation(s) with missing Block.IDs or Treatment were removed.")
    Block.IDs <- Block.IDs[-ind.na]
    Treatment <- Treatment[-ind.na]
    for (i in 1:length(Ks)) {
      Ks[[i]] <- Ks[[i]][-ind.na, -ind.na]
    }
  }

  ind <- order(Block.IDs, Treatment)
  Block.IDs <- Block.IDs[ind]
  Treatment <- Treatment[ind]
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][ind, ind]
  }

  treat.levels <- sort(unique(Treatment))
  lab <- match(Treatment, treat.levels) - 1L

  m <- length(treat.levels)
  N <- length(Treatment)

  if (m == 1) {
    stop("at least two treatment labels are needed.")
  }
  if (max(table(Block.IDs)) < 2) {
    stop("at least one block with multiple observations is needed.")
  }

  # Centering, symmetrization and normalization of the kernels (these do not change the P-values)

  for (i in 1:length(Ks)) {
    nK <- Ks[[i]] - (as.matrix(rep(1, N)) %*% t((1/N) * (Ks[[i]] %*% as.matrix(rep(1, N))))) - (((1/N) * (Ks[[i]] %*% as.matrix(rep(1, N)))) %*% rep(1, N)) + as.numeric((1/N^2) * (rep(1, N) %*% Ks[[i]] %*% as.matrix(rep(1, N))))
    nK <- (nK + t(nK))/2
    Ks[[i]] <- nK/sum(diag(nK))
  }
  names(Ks) <- kern.names

  # Test statistic values T = tr(X^T K X) (Eq. 7) for the original data (first row) and the rearranged
  # data (the other rows), for each kernel (columns), computed in C++ (src/kernrbd.cpp), where the
  # rearrangements are identical to those by base::sample for the same random seed

  T.null <- function(blocks, use) {
    kernrbd_null_T(Ks, lab, lapply(blocks, function(x) as.integer(x) - 1L), as.integer(use), n.res, n.cores)
  }

  # P-values of all elements of a reference set, where the first element is the original data:
  # P_(r) = (1/(R + 1)) * sum_{r' = 0}^{R} I(T_(r') >= T_(r)), r = 0, ..., R (Eqs. 8, 9, 14 and 15).
  # Test statistic values that differ only by floating-point errors are treated as ties.

  ref.pvals <- function(T.all) {
    R1 <- length(T.all)
    eps <- sqrt(.Machine$double.eps) * max(1, abs(T.all))
    (R1 - findInterval(T.all - eps, sort(T.all), left.open = TRUE))/R1
  }

  # Omnibus P-values of all elements of a reference set, given the minimum P-values:
  # P_omnibus(r) = (1/(R + 1)) * sum_{r' = 0}^{R} I(M_(r') <= M_(r)), r = 0, ..., R (Eqs. 16 and 17)

  ref.omni.pvals <- function(M.all) {
    rank(M.all, ties.method = "max")/length(M.all)
  }

  # Step-down minP adjusted P-values (Eq. 10), where P.all is the (R + 1) x q matrix of
  # P-values and its first row is for the original data

  stepdown.minP <- function(P.all) {
    q <- ncol(P.all)
    p.obs <- P.all[1, ]
    o <- order(p.obs)
    P.ord <- P.all[, o, drop = FALSE]
    Q <- P.ord
    if (q > 1) {
      for (l in (q - 1):1) {
        Q[, l] <- pmin(Q[, l + 1], P.ord[, l])
      }
    }
    adj <- sapply(1:q, function(l) mean(Q[, l] <= p.obs[o][l]))
    adj <- cummax(adj)
    adj[order(o)]
  }

  # Global testing: within-block rearrangement of all treatment labels, where the same
  # rearrangements are used for all kernels

  T.global <- T.null(split(1:N, Block.IDs), 0:(m - 1))
  global.P.all <- apply(T.global, 2, ref.pvals)
  global.P.all <- matrix(global.P.all, ncol = length(Ks))

  itembyitem.pvals <- global.P.all[1, ]
  names(itembyitem.pvals) <- kern.names
  minP.pval <- ref.omni.pvals(apply(global.P.all, 1, min))[1]
  global.out <- list(itembyitem.pvals = itembyitem.pvals, minP.pval = minP.pval)

  if (m == 2) {
    output <- list(Ks = Ks, block.IDs = Block.IDs, treatment = Treatment, global.out = global.out)
  }

  if (m > 2) {

    # Pairwise comparisons: within-block rearrangement of the two treatment labels only,
    # where the same rearrangements are used for all kernels

    pair.grps <- t(combn(treat.levels, 2))
    n.pairs <- nrow(pair.grps)
    block.levels <- unique(Block.IDs)

    pair.P.all <- lapply(1:length(Ks), function(i) matrix(NA, n.res + 1, n.pairs))
    for (j in 1:n.pairs) {
      treat.ind <- is.element(Treatment, pair.grps[j, ])
      treat.block.inds <- lapply(block.levels, function(b) which(treat.ind & Block.IDs == b))
      treat.block.inds <- treat.block.inds[sapply(treat.block.inds, length) > 1]
      T.pair <- T.null(treat.block.inds, match(pair.grps[j, ], treat.levels) - 1L)
      for (i in 1:length(Ks)) {
        pair.P.all[[i]][, j] <- ref.pvals(T.pair[, i])
      }
    }

    conc.pairs <- paste(pair.grps[, 1], pair.grps[, 2], sep = "-")

    itembyitem.pvals <- t(sapply(pair.P.all, function(P) P[1, ]))
    itembyitem.adj.pvals <- t(sapply(pair.P.all, stepdown.minP))
    itembyitem.pvals <- matrix(itembyitem.pvals, nrow = length(Ks), dimnames = list(kern.names, conc.pairs))
    itembyitem.adj.pvals <- matrix(itembyitem.adj.pvals, nrow = length(Ks), dimnames = list(kern.names, conc.pairs))

    # Omnibus test for each pair

    omni.P.all <- sapply(1:n.pairs, function(j) {
      ref.omni.pvals(apply(sapply(pair.P.all, function(P) P[, j]), 1, min))
    })
    omni.P.all <- matrix(omni.P.all, ncol = n.pairs)
    minP.pvals <- omni.P.all[1, ]
    minP.adj.pvals <- stepdown.minP(omni.P.all)
    names(minP.pvals) <- conc.pairs
    names(minP.adj.pvals) <- conc.pairs

    pairwise.out <- list(itembyitem.pvals = itembyitem.pvals, itembyitem.adj.pvals = itembyitem.adj.pvals,
                         minP.pvals = minP.pvals, minP.adj.pvals = minP.adj.pvals)

    output <- list(Ks = Ks, block.IDs = Block.IDs, treatment = Treatment, global.out = global.out, pairwise.out = pairwise.out)

  }

  attr(output, "n.res") <- n.res
  attr(output, "n.removed") <- length(ind.na)
  class(output) <- "KernRBD"
  return(output)
}
