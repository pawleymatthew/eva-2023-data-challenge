# update step
cp_decomp_kz_update <- function(S, i) {
  d <- ncol(S)
  if (d == 1) {
    S_update <- NULL
    cp_col <- sqrt(S)
    D_i <- NULL
  } else {
    D_i <- max((1 / S[i, i]) * (S[i, -i] %*% t(S[i, -i])) * (1 / S[-i, -i]))
    M_i <- max(D_i, 1)
    a <- sapply(1:d, function(j) {
      ifelse(i == j, 
             sqrt(S[i, i]) * M_i, 
             S[j, i] / sqrt(S[i, i] * M_i))
    })
    cp_col <- a
    S_update <- S[-i, -i] - (cp_col[-i] %*% t(cp_col[-i]))
  }
  return(list("cp_col" = cp_col, "S_update" = S_update, "D_i" = D_i))
}

# full algorithm for a given path of length d
cp_decomp_kz <- function(S, path = 1:ncol(S)) {
  d <- ncol(S)
  if (pracma::Rank(S) < d) stop("Matrix S is not full-rank.")
  # preallocate
  A <- matrix(0, nrow = d, ncol = d)
  rownames(A) <- rownames(S)
  D <- rep(NA, d - 1)
  # repeatedly apply update step
  S_current <- S
  i_current <- path[1]
  for (iter in 1:d) {
    temp <- cp_decomp_kz_update(S = S_current, i = i_current)
    A[sort(path[iter:d]), iter] <- temp$cp_col
    if (iter < d) D[iter] <- temp$D_i
    S_current <- temp$S_update
    i_current <- rank(path[(iter + 1):d])[1]
  }
  err <- norm((A %*% t(A)) - S, "F")
  if (is.nan(err)) {err <- Inf}
  return(list("A" = A, "D" = D, "err" = err))
}

# search the path space for a decomposition
# either sample paths randomly, or choose a start index, e.g. path_start = 2
cp_decomp_kz_search <- function(S, path_start = "random", max_paths = 10) {
  d <- ncol(S)
  cp_success <- FALSE
  iter <- 0
  start_time <- Sys.time()
  while(!cp_success & (iter < max_paths)) {
    if (is.numeric(path_start)) {
      if (path_start > d) stop("path_start must be between 1 and ncol(S)")
      path <- c(path_start, sample(c(1:d)[-path_start], d - 1))
    } else {
      path <- sample(1:d, d)
    } 
    temp <- cp_decomp_kz(S = S, path = path)
    cp_success <- (temp$err < 1e-12)
    iter <- iter + 1
  }
  elapsed_time <- Sys.time() - start_time
  if (cp_success) {
    return(list("A" = temp$A, "path" = path, "exact_decomp" = TRUE, "n_paths" = iter, "elapsed_time" = elapsed_time))
  } else {
    return(list("A" = NULL, "path" = NULL, "exact_decomp" = FALSE, "n_paths" = iter, "elapsed_time" = elapsed_time))
  }
}