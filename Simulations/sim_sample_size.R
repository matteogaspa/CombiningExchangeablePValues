rm(list = ls())
library(gtools)
source("functions_pvals.R")

# function to obtain p-values -----
simul_pvals <- function(n, n0, K, mu){
  sim_0 <- rnorm(n0, mean = mu)
  mu_0  <- (1/sqrt(n0)) * sum(sim_0)
  mu_i  <- 1:K
  t_i   <- 1:K
  p_i   <- 1:K
  for(i in 1:K){
    sim_i   <- rnorm(n[i], mean = mu)
    mu_i[i] <- (1/sqrt(n[i])) * sum(sim_i)
    t_i[i]  <- (mu_i[i] + mu_0)/sqrt(2)
    p_i[i]  <- 2 * min(pnorm(t_i[i]), pnorm(-t_i[i]))
  }
  return(list(
    "s_tets" = t_i,
    "p_vals" = p_i
  ))
}

# set the parameters -----
n0 <- 25
K  <- 10
n  <- 10 * 1:10
mu <- 0

# simulation under H0 -----
B <- 10000
mat_t <- mat_p <- matrix(NA, B, K)
mat_avg <- mat_med <- matrix(NA, B, 4)
set.seed(123)
for(i in 1:B){
  res       <- simul_pvals(n, n0, K, mu)
  mat_t[i,] <- res$s_tets
  mat_p[i,] <- res$p_vals
  
  # merging p-vals
  
  # decreasing sample size
  mat_avg[i,1] <- f_e_bavg(mat_p[i,])
  mat_med[i,1] <- f_e_ruger(mat_p[i,], K/2)
  
  # increasing sample size
  mat_avg[i,2] <- f_e_bavg(mat_p[i,K:1])
  mat_med[i,2] <- f_e_ruger(mat_p[i,K:1], K/2)
  
  # random order 
  mat_avg[i,3] <- f_e_bavg(permute(mat_p[i,]))
  mat_med[i,3] <- f_e_ruger(permute(mat_p[i,]), K/2)
  
  # standard rule
  mat_avg[i, 4] <- f_bavg(mat_p[i,])
  mat_med[i, 4] <- f_ruger(mat_p[i,], K/2)
}

cor(mat_t)

par(mfrow=c(2,5))
for(i in 1:K){
  plot.ecdf(mat_p[,i])
}

par(mfrow=c(1,3))
for(i in 1:3){
  plot.ecdf(mat_avg[,i])
  abline(0, 1, col = "red", lty = 2)
}
# same results under any order

for(i in 1:3){
  plot.ecdf(mat_med[,i])
  abline(0, 1, col = "red", lty = 2)
}
# same results under any order

# Simulation study ------
B     <- 10000                     # number of replications
mus   <- seq(0, 0.5, by = 0.01)    # tested means
alpha <- 0.05                      # confidence level

res_avg <- res_med <- matrix(NA, nrow = length(mus), ncol = 4)

set.seed(1234)
for(j in 1:length(mus)){
  mat_avg <- matrix(NA, nrow = B, ncol = 4)
  mat_med <- matrix(NA, nrow = B, ncol = 4)
  for(i in 1:B){
    res    <- simul_pvals(n, n0, K, mus[j])
    z_test <- res$s_tets
    p_vals <- res$p_vals
    
    # merging p-vals
    
    # decreasing sample size
    mat_avg[i,1] <- f_e_bavg(p_vals) 
    mat_med[i,1] <- f_e_ruger(p_vals, K/2)
    
    # increasing sample size
    mat_avg[i,2] <- f_e_bavg(p_vals[K:1])
    mat_med[i,2] <- f_e_ruger(p_vals[K:1], K/2)
    
    # random order 
    mat_avg[i,3] <- f_e_bavg(permute(p_vals))
    mat_med[i,3] <- f_e_ruger(permute(p_vals), K/2)
    
    # standard rules
    mat_avg[i,4] <- f_bavg(p_vals)
    mat_med[i,4] <- f_ruger(p_vals, K/2)
  }
  
  res_avg[j,] <- colMeans(I(mat_avg <= alpha))
  res_med[j,] <- colMeans(I(mat_med <= alpha))
  
  cat("Iter", j, "\n")
}

# plots -----
#par(las = 1, family = "serif", font.main = 1.5)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))

par(mar = c(4, 4, 2, 1))
plot(mus, smooth.spline(mus, res_avg[,1])$y, type = "l", ylab = "power",
     xlab = expression(mu), ylim = c(0,1), lwd = 3, lty = 2,
     main = " 'Twice the average' ", cex.lab = 1.5, cex.main = 1.5)
lines(mus, smooth.spline(mus, res_avg[,2])$y, col = "red", lty = 1, lwd = 2)
lines(mus, smooth.spline(mus, res_avg[,3])$y, col = "blue", lty = 3, lwd = 3)
lines(mus, smooth.spline(mus, res_avg[,4])$y, col = "orange", lty = 4, lwd = 3)
#legend("bottom", inset = 0.01, c("Increasing order", "Decreasing order", "Random order", "Usual"), horiz = T, text.width = strwidth("Q"),
#       lty = c(2, 1, 3, 4), lwd = c(1.5, 1, 3, 1.5), col = c("black", "red", "blue", "orange"), box.col = "white", cex = 0.8, bg = "white")
abline(h=0.05, col = "gray", lty = 5)
text(x = 0.45, y = 0.1, labels = expression(alpha == 0.05), col = "gray", cex = 1)

par(mar = c(4, 4, 2, 1))
plot(mus, smooth.spline(mus, res_med[,1])$y, type = "l", ylab = "power",
     xlab = expression(mu), ylim = c(0,1), lwd = 3, lty = 2,
     main = " 'Twice the median' ", cex.lab = 1.5, cex.main = 1.5)
lines(mus, smooth.spline(mus, res_med[,2])$y, col = "red", lty = 1, lwd = 2)
lines(mus, smooth.spline(mus, res_med[,3])$y, col = "blue", lty = 3, lwd = 3)
lines(mus, smooth.spline(mus, res_med[,4])$y, col = "orange", lty = 4, lwd = 3)
#legend("bottom", inset = 0.01, c("Increasing order", "Decreasing order", "Random order", "Usual"), horiz = T, text.width = strwidth("Q"),
#       lty = c(2, 1, 3, 4), lwd = c(1.5, 1, 3, 1.5), col = c("black", "red", "blue", "orange"), box.col = "white", cex = 0.8, bg = "white")
abline(h=0.05, col = "gray", lty = 5)
text(x = 0.45, y = 0.1, labels = expression(alpha == 0.05), col = "gray", cex = 1)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", c("Increasing order", "Decreasing order", "Random order", "Non ex-p"), x.intersp = 0.15, text.width = 0.185,
       lty = c(2, 1, 3, 4), lwd = c(3, 2, 3, 3), col = c("black", "red", "blue", "orange"), box.col = "gray", cex = 1.25, bg = "white", horiz = T)

