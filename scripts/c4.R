rm(list = ls())

library(cluster)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)
library(evd)
library(dplyr)
library(purrr)
library(lattice)
library(gridExtra)

# Directory -------------------------------------------------------------------------

data_dir <- file.path("data")
fig_dir <- file.path("figs", "c4")
fn_dir <- file.path("R")

# Source functions and load data ----------------------------------------------------

fn_files <- list.files(path = fn_dir, pattern = "*.R$", recursive = TRUE, full.names = TRUE)
sapply(fn_files, source)

# Load data -------------------------------------------------------------------------

UtopulaU1 <- read.csv(file = file.path(data_dir, "UtopulaU1.csv"), header = TRUE)
UtopulaU2 <- read.csv(file = file.path(data_dir, "UtopulaU2.csv"), header = TRUE)

colnames(UtopulaU1) <- paste(colnames(UtopulaU1), "U1", sep = ".")
colnames(UtopulaU2) <- paste(colnames(UtopulaU2), "U2", sep = ".")

Utopula <- cbind(UtopulaU1, UtopulaU2)
Utopula <- as.matrix(Utopula)

n <- nrow(Utopula)
d <- ncol(Utopula)

Utopula <- exp(Utopula / 2) # Frechet(2) margins

# C4 --------------------------------------------------------------------------------

k <- 250
K <- 5

u1 <- function(varnames) dplyr::case_when(endsWith(varnames, "U1") ~ evd::qfrechet(299/300, shape = 2),
                                          endsWith(varnames, "U2") ~ evd::qfrechet(288/300, shape = 2))
u2 <- evd::qfrechet(299/300, shape = 2)

## Clustering ------------------------------------------------------------------------

F_Utopula <- evd::pfrechet(Utopula, shape = 2)
D <- dist(t(F_Utopula), method = "manhattan") / (2 * n)
Utopula_cluster <- cluster::pam(x = D, k = K, diss = TRUE, cluster.only = TRUE)
cluster_dims <- table(Utopula_cluster)

X_clusters <- lapply(1:K, function(l) Utopula[, Utopula_cluster == l])

## Empirical estimation approach -----------------------------------------------------

# empirical A for each cluster
A_hat_clusters <- lapply(X_clusters, function(X_l) max_linear_A_emp(X_l, k, alpha = 2))

p_u1_emp_clusters <- lapply(A_hat_clusters, function(A) max_linear_p_hat(A = A, u = u1(rownames(A)), alpha = 2)) 
p_u1_emp <- p_u1_emp_clusters %>% purrr::reduce(prod)

p_u2_emp_clusters <- lapply(A_hat_clusters, function(A) max_linear_p_hat(A = A, u = u2, alpha = 2)) 
p_u2_emp <- p_u2_emp_clusters %>% purrr::reduce(prod)

# Sparse empirical estimation approach ----------------------------------------------

Utopula1 <- Utopula^2 # Frechet shape 1

X1_clusters <- lapply(1:K, function(l) Utopula1[, Utopula_cluster == l])
A_hatstar_clusters <- lapply(X1_clusters, function(X_l) {
  R_l <- apply(X_l, 1, pracma::Norm, p = 1)
  R_l_kplus1 <- Rfast::nth(R_l, k = k + 1, descending = TRUE)
  ext_ind <- which(R_l > R_l_kplus1)
  piTheta_ext_l <- (X_l[ext_ind, ] / R_l_kplus1) %>% apply(1, euc_proj) %>% t()
  d_l <- ncol(X_l)
  A_hatstar_l <- (d_l / k) *  t(piTheta_ext_l)
  return(A_hatstar_l)
})

p_u1_spemp_clusters <- lapply(A_hatstar_clusters, function(A) max_linear_p_hat(A = A, u = u1(rownames(A))^2, alpha = 1)) 
p_u1_spemp <- p_u1_spemp_clusters %>% purrr::reduce(prod)

p_u2_spemp_clusters <- lapply(A_hatstar_clusters, function(A) max_linear_p_hat(A = A, u = u2^2, alpha = 1)) 
p_u2_spemp <- p_u2_spemp_clusters %>% purrr::reduce(prod)

## CP estimation approach ------------------------------------------------------------

# estimate TPDM for each cluster
Sigma_hat_clusters <- lapply(X_clusters, function(X_l) tpdm_estimate(X_l, k))
# summary stats of TPDMs
Sigma_hat_clusters %>% lapply(FUN = function(S) unlist(S[upper.tri(S, diag = F)])) %>% lapply(summary) 

Theta_hat_kz <- tibble("cluster" = rep(1:K, cluster_dims), "path_start" = sequence(cluster_dims)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate("cp_decomp_kz" = list(cp_decomp_kz_search(Sigma_hat_clusters[[cluster]], path_start = path_start, max_paths = 100))) %>%
  tidyr::unnest_wider(cp_decomp_kz) %>%
  dplyr::rowwise() %>%
  dplyr::mutate("p_u1" = max_linear_p_hat(A = A, beta = 1:nrow(A), u = u1(names(which(Utopula_cluster == cluster))), alpha = 2)) %>%
  dplyr::mutate("p_u2" = max_linear_p_hat(A = A, u = u2, alpha = 2)) %>%
  dplyr::ungroup()

p_u1_cp <- Theta_hat_kz %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::reframe("p_u1" = list(unique(p_u1))) %>% 
  dplyr::pull(p_u1) %>% 
  expand.grid() %>% 
  apply(1, prod)

p_u2_cp <- Theta_hat_kz %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::reframe("p_u2" = list(unique(p_u2))) %>% 
  dplyr::pull(p_u2) %>% 
  expand.grid() %>% 
  apply(1, prod)

# Plots -----------------------------------------------------------------------------

colpal <- colorRampPalette(c('white', 'black'))(30)
breaks <- lattice::do.breaks(c(0, 1.3), length(colpal))

fig_name <- "Sigma-hat-12-levelplot.pdf"
pdf(file.path(fig_dir, fig_name), width = 10, height = 5)
S1 <- t(Sigma_hat_clusters[[1]][nrow(Sigma_hat_clusters[[1]]):1, ])
rownames(S1) <- paste0("X", 1:50)[Utopula_cluster == 1]
colnames(S1) <- rev(rownames(S1))
p1 <- lattice::levelplot(S1, col.regions = colpal, at = breaks, 
                   xlab = "", ylab = "",
                   colorkey = TRUE,
                   scales = list(x = list(cex = 1.1), y = list(cex = 1.1)))
S2 <- t(Sigma_hat_clusters[[2]][nrow(Sigma_hat_clusters[[2]]):1, ])
rownames(S2) <- paste0("X", 1:50)[Utopula_cluster == 2]
colnames(S2) <- rev(rownames(S2))
p2 <- lattice::levelplot(S2, col.regions = colpal, at = breaks, 
                   xlab = "", ylab = "",
                   colorkey = TRUE,
                   scales = list(x = list(cex = 1.1), y = list(cex = 1.1)))
grid.arrange(p1, p2, ncol = 2)
dev.off()

fig_name <- "prob-estimates.pdf"
pdf(file.path(fig_dir, fig_name), width = 10, height = 6)
par(mfrow = c(1,2))
#left
y <- log10(c(unlist(p_u1_emp_clusters), unlist(p_u1_spemp_clusters), Theta_hat_kz$p_u1))
plot(c(1:5, 1:5, Theta_hat_kz$cluster), y, type = "n", xlab = "", ylab = "")
abline(v = 1:5, lty = 2, col = "grey")
points(Theta_hat_kz$cluster, log10(Theta_hat_kz$p_u1), pch = 4, cex = 1.5)
points(1:5, log10(unlist(p_u1_emp_clusters)), pch = 20, cex = 1.5)
points(1:5, log10(unlist(p_u1_spemp_clusters)), pch = 17, cex = 1.5)
legend("topright",
       legend = c("CP", "Empirical", "Sparse emp.     "), pch = c(4, 20, 17))
title(xlab = "Cluster", ylab = expression(log[10] ~ probability))
# right
boxwex <- 0.4
boxplot(log10(p_u1_cp), 
        log10(p_u2_cp), 
        names = c(expression(bold(u) == bold(u)[1]), expression(bold(u) == bold(u)[2])),
        xlab = "Threshold vector",
        ylab = expression(log[10] ~ probability),
        boxwex = boxwex)
dev.off()

cat("Empirical estimate of p1:", p_u1_emp, "\n")
cat("Empirical estimate of p2:", p_u2_emp, "\n")
cat("Sparse empirical estimate of p1:", p_u1_spemp, "\n")
cat("Sparse empirical estimate of p2:", p_u2_spemp, "\n")
cat("CP estimate of p1:", median(p_u1_cp), "\n")
cat("CP estimate of p2:", median(p_u2_cp), "\n")


