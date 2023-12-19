rm(list = ls())

library(dplyr)
library(lattice)
library(gridExtra)
library(Ternary)
library(Rfast)
library(stringr)

# Directory -------------------------------------------------------------------------

data_dir <- file.path("data")
fig_dir <- file.path("figs", "c3")
fn_dir <- file.path("R")

# Source functions and load data ----------------------------------------------------

fn_files <- list.files(path = fn_dir, pattern = "*.R$", recursive = TRUE, full.names = TRUE)
sapply(fn_files, source)

# Load data -------------------------------------------------------------------------

Coputopia <- read.csv(file = file.path(data_dir, "Coputopia.csv"), header = TRUE)

X <- Coputopia %>% dplyr::select(Y1, Y2, Y3) %>% as.matrix() %>% exp() # Frechet(1) margins
colnames(X) <- paste0("X", seq_len(ncol(X)))
d <- ncol(X)

# C3 --------------------------------------------------------------------------------

k <- 500 # value used for submission

R <- apply(X, 1, pracma::Norm, p = 1)
R_kplus1 <- Rfast::nth(R, k = k + 1, descending = TRUE)
ext_ind <- which(R > R_kplus1)
piTheta_ext <- (X[ext_ind, ] / R_kplus1) %>% apply(1, euc_proj) %>% t()

face <- apply(piTheta_ext, 1, function(x) stringr::str_c(as.character(which(x > 0)), collapse = ""))
face_counts <- table(face)

A_hat <- max_linear_A_emp(X, k, alpha = 1)
A_hat_star <- (d / k) *  t(piTheta_ext)


## C3.1 ------------------------------------------------------------------------------

beta1 <- c(1, 2, 3)
u1 <- exp(6)
p1_hat <- max_linear_p_hat(A_hat_star, beta1, u1, alpha = 1)

## C3.2 ------------------------------------------------------------------------------

beta2 <- c(1, 2)
u2 <- exp(7)
p2_hat <- max_linear_p_hat(A_hat_star, beta2, u2, alpha = 1)

# Sensitivity to choice of k  --------------------------------------------------------

k_vals <- seq(from = 50, to = 1500, by = 25)
pk <- array(NA, dim = c(length(k_vals), 3))
colnames(pk) <- c("k", "p1", "p2")
pk[, "k"] <- k_vals

for (i in seq_along(k_vals)) {
  R_kplus1_k <- Rfast::nth(R, k = k_vals[i] + 1, descending = TRUE)
  ext_ind_k <- which(R > R_kplus1_k)
  piTheta_ext_k <- (X[ext_ind_k, ] / R_kplus1_k) %>% apply(1, euc_proj) %>% t()
  Ak <- (d / k_vals[i]) *  t(piTheta_ext_k)
  pk[i, "p1"] <- max_linear_p_hat(Ak, beta1, u1, alpha = 1)
  pk[i, "p2"] <- max_linear_p_hat(Ak, beta2, u2, alpha = 1)
}
pk <- as.data.frame(pk)

# Figures ---------------------------------------------------------------------------

colpal <- colorRampPalette(c('white', 'black'))(30)
ncolA <- 100 # number of columns of the A estimates to plot

# levelplot of A_hat
fig_name <- "A.pdf"
pdf(file.path(fig_dir, fig_name), width = 10, height = 3)
par(mar = c(0, 0, 0, 0))
breaks <- lattice::do.breaks(c(0, 0.006), length(colpal))
# left
p1 <- levelplot(t(A_hat[3:1, 1:ncolA]), col.regions = colpal, at = breaks, 
          main = "", xlab = "", ylab = "", 
          colorkey = list(labels = list(cex = 1.2)),
          scales = list(x = list(at = c(1,100), labels = c(1,100), cex = 1.4),
                        y = list(labels = c("3", "2", "1"), cex = 1.4))) %>% 
  update(aspect = 0.5)
# right
p2 <- levelplot(t(A_hat_star[3:1, 1:ncolA]),  col.regions = colpal, at = breaks, 
                main = "", xlab = "", ylab = "",
                colorkey = list(labels = list(cex = 1.2)),
                scales = list(x = list(at = c(1,100), labels = c(1,100), cex = 1.4),
                              y = list(labels = c("3", "2", "1"), cex = 1.4))) %>% 
  update(aspect = 0.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# ternary plots
fig_name <- "ternary.pdf"
pdf(file.path(fig_dir, fig_name), width = 10, height = 5)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 2))
# left
TernaryPlot(atip = "{1}", btip = "{2}", ctip = "{3}",
            alab = "{1,3}", blab = "{1,2}", clab = "{2,3}", lab.font = 2, lab.offset = 0.1,
            lab.cex = 1.4, tip.cex = 1.4,
            grid.lines = 0)
TernaryPoints((k / d) * t(A_hat), col = "black", pch = 4)
TernaryText(c(1, 1, 1), labels = c("{1,2,3}"), col = "black", font = 2, cex = 1.4)
# right
TernaryPlot(atip = "{1}", btip = "{2}", ctip = "{3}",
            alab = "{1,3}", blab = "{1,2}", clab = "{2,3}", lab.font = 2, lab.offset = 0.1,
            lab.cex = 1.4, tip.cex = 1.4,
            grid.lines = 0)
TernaryPoints(piTheta_ext[face == "123", ], col = "black", pch = 4)
TernaryPoints(piTheta_ext[face %in% c("12", "13", "23"), ], col = "red", pch = 4)
TernaryPoints(piTheta_ext[face %in% c("1", "2", "3"), ], col = "blue", pch = 4)
TernaryText(c(1, 1, 1), labels = c("{1,2,3}"), col = "black", font = 2, cex = 1.4)
dev.off()

fig_name <- "k-sensitivity.pdf"
pdf(file.path(fig_dir, fig_name), width = 8, height = 5)
par(mfrow = c(1, 1))
plot(k_vals, log10(pk[, "p1"]), col = "red", type = "l", lwd = 1.2,
     xlab = expression(k), 
     ylab = expression(log[10] ~ probability),
     xlim = range(k_vals))
abline(v = 500, lty = 2) # my k
abline(h = log10(5.38e-5), lty = 2, col = "red") # true value
points(k_vals, log10(pk[, "p2"]), col = "blue", type = "l", lwd = 1.2)
abline(h = log10(2.98e-5), lty = 2, col = "blue") # true value
legend("bottomright", 
       legend = c(expression(Estimate ~ hat(p)[1]), 
                  expression(Estimate ~ hat(p)[2]), 
                  expression(True ~ p[1]), 
                  expression(True ~ p[2])),
       lty = c(1, 1, 2, 2), col = c("red", "blue", "red", "blue"))
dev.off()

