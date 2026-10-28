D.to.K <-
function(D) {
  if (is.list(D) && !inherits(D, "dist")) {
    return(lapply(D, D.to.K))
  }
  D <- as.matrix(D)
  if (!is.numeric(D) || nrow(D) != ncol(D)) {
    stop("D should be a numeric square distance matrix (or a dist object, or a list of them).")
  }
  n <- nrow(D)
  centerM <- diag(n) - 1/n
  K <- -0.5 * centerM %*% (D * D) %*% centerM
  dimnames(K) <- dimnames(D)
  return(K)
}
