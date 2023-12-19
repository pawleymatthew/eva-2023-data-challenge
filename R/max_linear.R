# \hat{A}_k, the empirical estimate of the noise coefficient matrix
max_linear_A_emp <- function(X, k, alpha) {
  R <- apply(X, 1, pracma::Norm, p = alpha)
  R_k <- Rfast::nth(R, k = k, descending = TRUE)
  ext_ind <- R >= R_k
  Theta_ext <- X[ext_ind, ] / R[ext_ind]
  A_hat <- (ncol(X) / k)^(1 / alpha) * t(Theta_ext)
  return(A_hat)
}

# \hat{\mathbb{P}}(X \in C_{beta, u})
max_linear_p_hat <- function(A, beta = 1:nrow(A), u, alpha = 2) {
  beta_cols <- apply(A, 2, function(x) (all(x[beta] > 0) & all(x[-beta] == 0)))
  p <- apply(A[, beta_cols, drop = FALSE], 2, function(a) min(a[beta] / u)^alpha) %>% sum()
  return(p)
}

# empirical TPDM, \hat{\Sigma}
tpdm_estimate <- function(X, k) {
  A <- max_linear_A_emp(X, k, alpha = 2)
  Sigma_hat <- A %*% t(A)
  return(Sigma_hat)
}

# Euclidean projection onto L1-simplex
euc_proj <- function(x) {
  b <- 1 # radius of the simplex
  u <- sort(x, decreasing = TRUE)
  sx <- cumsum(u)
  rho <- which(u > (sx - b) / (1:length(x)))
  theta <- max(0, (sx[rho] - b) / rho)
  w <- x - theta
  w[w < 0] = 0
  return(w)
}