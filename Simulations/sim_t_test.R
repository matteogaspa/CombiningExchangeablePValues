rm(list = ls())
source(file = "functions_pvals.R")

# Functions
generate_sample <- function(n, mus, sigmas){
  # mus is the mean and sigmas is the vector containing the variance of the rvs
  # return a matrix n x K
  K    <- length(mus)
  samp <- matrix(NA, nrow = n, ncol = K)
  for(i in 1:K){
    samp[,i] <- rnorm(n, mean = mus[i], sd = sigmas[i])
  }
  return(samp)
}

test_null <- function(data){
  n         <- nrow(data)
  hat_mu    <- colMeans(data)
  hat_sigma <- apply(data, 2, var)
  t_tests   <- sqrt(n) * (hat_mu/sqrt(hat_sigma))
  p_vals    <- 2 * pt(-abs(t_tests), df = n-1)
  s_k       <- apply(data, 2, function(x) sum(x^2)/n)
  return(list(
    p_vals = p_vals,
    t_k    = t_tests,
    s_k    = s_k
  ))
}

# example under the null -----
B      <- 10000
n      <- 10
K      <- 20
mus    <- rep(0, K)
sigmas <- rep(1, K)

p_vals <- s_k <- t_k <- matrix(NA, nrow = B, ncol = K)
mean_p <- matrix(NA, nrow = B, ncol = 3)
medi_p <- matrix(NA, nrow = B, ncol = 3)
bonf_p <- matrix(NA, nrow = B, ncol = 2)
set.seed(10)
for(i in 1:B){
  camp       <- generate_sample(n, mus, sigmas)
  test_t     <- test_null(camp)
  p_vals[i,] <- test_t$p_vals
  t_k[i,]    <- test_t$t_k
  s_k[i,]    <- test_t$s_k
  p          <- test_t$p_vals
  # combine p-vals
  mean_p[i,1] <- f_e_bavg(p[order(test_t$s_k, decreasing = T)])
  mean_p[i,2] <- f_e_bavg(p[order(test_t$s_k, decreasing = F)])
  mean_p[i,3] <- f_e_bavg(p[sample(1:K)])
  
  medi_p[i,1] <- f_e_ruger(p[order(test_t$s_k, decreasing = T)], K/2)
  medi_p[i,2] <- f_e_ruger(p[order(test_t$s_k, decreasing = F)], K/2)
  medi_p[i,3] <- f_e_ruger(p[sample(1:K)], K/2)
  
  bonf_p[i,1] <- K*min(p)
  bonf_p[i,2] <- pchisq(-2*log(prod(p)), df = 2*K)
}

par(mfrow = c(2, 5))
for(i in 1:K){
  plot.ecdf(p_vals[,i])
  abline(0, 1, col = "gray", lty = 2)
}

cor(p_vals)

par(mfrow = c(4, 2))
for(i in 1:3){
  plot.ecdf(mean_p[,i], xlim = c(0,1))
  abline(0, 1, col = "gray", lty = 2)
  plot.ecdf(medi_p[,i], xlim = c(0,1))
  abline(0, 1, col = "gray", lty = 2)
}
plot.ecdf(bonf_p[,1], xlim = c(0,1))
abline(0, 1, col = "gray", lty = 2)
plot.ecdf(bonf_p[,2], xlim = c(0,1))
abline(0, 1, col = "gray", lty = 2)

colMeans(mean_p<=0.05)
colMeans(medi_p<=0.05)
colMeans(bonf_p<= 0.05)

for(i in 1:10){
  cat(cor(p_vals[,i], s_k[,i]), cor(s_k[,i], abs(t_k[,i])),"\n")
}

# simulation 1: \mu_k=\mu*k ------
mus_sq <- seq(0, 0.1, length = 50)
B      <- 10000
n      <- 10
K      <- 20
sigmas <- rep(1, K)

res_mean <- matrix(NA, nrow = length(mus_sq), ncol = 4)
res_medi <- matrix(NA, nrow = length(mus_sq), ncol = 4)
res_bonf <- matrix(NA, nrow = length(mus_sq), ncol = 2)

set.seed(10)
for(j in 1:length(mus_sq)){
  mus    <- mus_sq[j]*(1:K)
  mean_p <- matrix(NA, nrow = B, ncol = 4)
  medi_p <- matrix(NA, nrow = B, ncol = 4)
  bonf_p <- matrix(NA, nrow = B, ncol = 2)
  for(i in 1:B){
    camp       <- generate_sample(n, mus, sigmas)
    test_t     <- test_null(camp)
    p          <- test_t$p_vals
    # combine p-vals
    mean_p[i,1] <- f_e_bavg(p[order(test_t$s_k, decreasing = T)])
    mean_p[i,2] <- f_e_bavg(p[order(test_t$s_k, decreasing = F)])
    mean_p[i,3] <- f_e_bavg(p[sample(1:K)])
    mean_p[i,4] <- f_bavg(p)
    
    medi_p[i,1] <- f_e_ruger(p[order(test_t$s_k, decreasing = T)], K/2)
    medi_p[i,2] <- f_e_ruger(p[order(test_t$s_k, decreasing = F)], K/2)
    medi_p[i,3] <- f_e_ruger(p[sample(1:K)], K/2)
    medi_p[i,4] <- f_ruger(p, K/2)
    
    bonf_p[i,1] <- K*min(p)
    bonf_p[i,2] <- pchisq(-2*log(prod(p)), df = 2*K, lower.tail = F)
  }
  res_mean[j,] <- apply(mean_p, 2, function(x) mean(I(x <= 0.05)))
  res_medi[j,] <- apply(medi_p, 2, function(x) mean(I(x <= 0.05)))
  res_bonf[j,] <- apply(bonf_p, 2, function(x) mean(I(x <= 0.05))) 
  cat("Iter:",j,"\n")
}

save(res_mean, res_medi, res_bonf, file = "res_simul_t0.RData")


# plots -----
rm(list = ls())
load("res_simul_t0.RData")
mus <- seq(0, 0.20, length = 50)

layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.2))
par(mar = c(4, 4, 2, 1))


plot(mus, smooth.spline(mus, res_mean[,1])$y, type = "l", ylab = "power",
     xlab = expression(mu), ylim = c(0,1), lwd = 3, lty = 2,
     main = " 'Twice the average' ")
lines(mus, smooth.spline(mus, res_mean[,2])$y, col = "red", lty = 1, lwd = 2)
lines(mus, smooth.spline(mus, res_mean[,3])$y, col = "blue", lty = 3, lwd = 3)
lines(mus, smooth.spline(mus, res_mean[,4])$y, col = "orange", lty = 4, lwd = 3)
lines(mus, smooth.spline(mus, res_bonf[,1])$y, col = "lightpink2", lty = 5, lwd = 3)
lines(mus, smooth.spline(mus, res_bonf[,2])$y, col = "forestgreen", lty = 6, lwd = 3)
#legend("bottom", inset = 0.01, c("Increasing order", "Decreasing order", "Random order", "Usual", "Bonferroni", "Fisher"), horiz = T, text.width = strwidth("Q"),
#       lty = c(2, 1, 3, 4), lwd = c(1.5, 1, 3, 1.5, 1.5, 1.5), col = c("black", "red", "blue", "orange", "purple", "forestgreen"), box.col = "white", cex = 0.8, bg = "white")
abline(h=0.05, col = "gray", lty = 5)
#text(x = 0.45, y = 0.075, labels = expression(alpha == 0.05), col = "gray", cex = 1)

plot(mus, smooth.spline(mus, res_medi[,1])$y, type = "l", ylab = "power",
     xlab = expression(mu), ylim = c(0,1), lwd = 3, lty = 2,
     main = " 'Twice the median' ")
lines(mus, smooth.spline(mus, res_medi[,2])$y, col = "red", lty = 1, lwd = 3)
lines(mus, smooth.spline(mus, res_medi[,3])$y, col = "blue", lty = 3, lwd = 3)
lines(mus, smooth.spline(mus, res_medi[,4])$y, col = "orange", lty = 4, lwd = 3)
lines(mus, smooth.spline(mus, res_bonf[,1])$y, col = "lightpink2", lty = 5, lwd = 3)
lines(mus, smooth.spline(mus, res_bonf[,2])$y, col = "forestgreen", lty = 6, lwd = 3)
#legend("bottom", inset = 0.01, c("Increasing order", "Decreasing order", "Random order", "Usual", "Bonferroni", "Fisher"), horiz = T, text.width = strwidth("Q"),
#       lty = c(2, 1, 3, 4), lwd = c(1.5, 1, 3, 1.5, 1.5, 1.5), col = c("black", "red", "blue", "orange", "purple", "forestgreen"), box.col = "white", cex = 0.8, bg = "white")
abline(h=0.05, col = "gray", lty = 5)

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", c("Decreasing order", "Non ex-p", "Increasing order", "Bonferroni", "Random order", "Fisher"), lwd = c(2, 2, 1.5, 2, 2, 2), lty = c(2, 4, 1, 5, 3, 6), col = c("black", "orange", "red", "lightpink2", "blue", "forestgreen"), box.col = "gray", box.lty = 1, cex = 1, bg = "white", ncol = 3)


