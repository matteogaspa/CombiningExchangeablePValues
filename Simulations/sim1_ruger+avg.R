rm(list = ls())
source("functions_pvals.R")

# Function to simulate p-values ------
sim_exc_pvals <- function(n, K, rho=0, mu=0){
  # n: number of reps
  # K: number of p-values
  # rho: corr. parm
  # mu: mean
  z   <- rnorm(1*n)
  z_k <- x_k <- matrix(rnorm(K*n), nrow = n, ncol = K)
  for(i in 1:n){
    x_k[i,] <- pnorm(rho*z[i] + sqrt(1-rho^2)*z_k[i,] - mu)
  }
  return(x_k)
}

# Set hyperparameters -----
mus   <- seq(0, 3, by = 0.1)  # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.9                  # corr

## Ruger combination rule -----
mat_rug <- matrix(NA, nrow = length(mus), ncol = 3)

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute ruger, e-ruger, u-ruger
  rug_cr   <- apply(mat_pvals, 1, function(x) f_ruger(x, k = K/2))
  rug_cr   <- sapply(rug_cr, function(x) min(x, 1))
  rug_e_cr <- apply(mat_pvals, 1, function(x) f_e_ruger(x, k = K/2))
  rug_e_cr <- sapply(rug_e_cr, function(x) min(x, 1))
  rug_u_cr <- rep(NA, B)
  for(i in 1:B){
    rug_u_cr[i] <- f_u_ruger(mat_pvals[i,], k = K/2, u = runif(1))
  }
  rug_u_cr <- sapply(rug_u_cr, function(x) min(x, 1))
  
  # power of a test alpha
  mat_rug[j,1] <- mean(I(rug_cr <= alpha))
  mat_rug[j,2] <- mean(I(rug_e_cr <= alpha))
  mat_rug[j,3] <- mean(I(rug_u_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_rug, file="mat_ruger.RData")


## Twice the average combination rule -----
mat_avg <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute avg
  avg_cr   <- apply(mat_pvals, 1, function(x) f_avg(x))
  avg_cr   <- sapply(avg_cr, function(x) min(x, 1))
  bavg_cr   <- apply(mat_pvals, 1, function(x) f_bavg(x))
  bavg_cr   <- sapply(bavg_cr, function(x) min(x, 1))
  
  avg_e_cr <- apply(mat_pvals, 1, function(x) f_e_avg(x))
  avg_e_cr <- sapply(avg_e_cr, function(x) min(x, 1))
  bavg_e_cr <- apply(mat_pvals, 1, function(x) f_e_bavg(x))
  bavg_e_cr <- sapply(bavg_e_cr, function(x) min(x, 1))
  
  avg_u_cr <- rep(NA, B)
  bavg_u_cr <- rep(NA, B)
  for(i in 1:B){
    u <- runif(1)
    avg_u_cr[i] <- f_u_avg(mat_pvals[i,], u = u)
    bavg_u_cr[i] <- f_u_bavg(mat_pvals[i,], u = u)
  }
  avg_u_cr <- sapply(avg_u_cr, function(x) min(x, 1))
  bavg_u_cr <- sapply(bavg_u_cr, function(x) min(x, 1))

  
  # power of a test alpha
  mat_avg[j,1] <- mean(I(avg_cr <= alpha))
  mat_avg[j,2] <- mean(I(avg_e_cr <= alpha))
  mat_avg[j,3] <- mean(I(avg_u_cr <= alpha))
  mat_avg[j,4] <- mean(I(bavg_cr <= alpha))
  mat_avg[j,5] <- mean(I(bavg_e_cr <= alpha))
  mat_avg[j,6] <- mean(I(bavg_u_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_avg,file= "mat_avg.RData")

## Bonferroni -----
mat_bonf <- rep(NA, length(mus))

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute ruger, e-ruger, u-ruger
  rug_cr   <- apply(mat_pvals, 1, function(x) f_ruger(x, k = 1))
  rug_cr   <- sapply(rug_cr, function(x) min(x, 1))
  
  # power of a test alpha
  mat_bonf[j] <- mean(I(rug_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_bonf, file = "mat_bonf.Rdata")

# Set hyperparameters (different rho) -----
mus   <- seq(0, 3, by = 0.1)  # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.1                  # corr


## Ruger combination rule -----
mat_rug <- matrix(NA, nrow = length(mus), ncol = 3)

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute ruger, e-ruger, u-ruger
  rug_cr   <- apply(mat_pvals, 1, function(x) f_ruger(x, k = K/2))
  rug_cr   <- sapply(rug_cr, function(x) min(x, 1))
  rug_e_cr <- apply(mat_pvals, 1, function(x) f_e_ruger(x, k = K/2))
  rug_e_cr <- sapply(rug_e_cr, function(x) min(x, 1))
  rug_u_cr <- rep(NA, B)
  for(i in 1:B){
    rug_u_cr[i] <- f_u_ruger(mat_pvals[i,], k = K/2, u = runif(1))
  }
  rug_u_cr <- sapply(rug_u_cr, function(x) min(x, 1))
  
  # power of a test alpha
  mat_rug[j,1] <- mean(I(rug_cr <= alpha))
  mat_rug[j,2] <- mean(I(rug_e_cr <= alpha))
  mat_rug[j,3] <- mean(I(rug_u_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_rug, file="mat_ruger2.RData")


## Twice the average combination rule -----
mat_avg <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute avg
  avg_cr   <- apply(mat_pvals, 1, function(x) f_avg(x))
  avg_cr   <- sapply(avg_cr, function(x) min(x, 1))
  bavg_cr   <- apply(mat_pvals, 1, function(x) f_bavg(x))
  bavg_cr   <- sapply(bavg_cr, function(x) min(x, 1))
  
  avg_e_cr <- apply(mat_pvals, 1, function(x) f_e_avg(x))
  avg_e_cr <- sapply(avg_e_cr, function(x) min(x, 1))
  bavg_e_cr <- apply(mat_pvals, 1, function(x) f_e_bavg(x))
  bavg_e_cr <- sapply(bavg_e_cr, function(x) min(x, 1))
  
  avg_u_cr <- rep(NA, B)
  bavg_u_cr <- rep(NA, B)
  for(i in 1:B){
    u <- runif(1)
    avg_u_cr[i] <- f_u_avg(mat_pvals[i,], u = u)
    bavg_u_cr[i] <- f_u_bavg(mat_pvals[i,], u = u)
  }
  avg_u_cr <- sapply(avg_u_cr, function(x) min(x, 1))
  bavg_u_cr <- sapply(bavg_u_cr, function(x) min(x, 1))
  
  
  # power of a test alpha
  mat_avg[j,1] <- mean(I(avg_cr <= alpha))
  mat_avg[j,2] <- mean(I(avg_e_cr <= alpha))
  mat_avg[j,3] <- mean(I(avg_u_cr <= alpha))
  mat_avg[j,4] <- mean(I(bavg_cr <= alpha))
  mat_avg[j,5] <- mean(I(bavg_e_cr <= alpha))
  mat_avg[j,6] <- mean(I(bavg_u_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_avg,file= "mat_avg2.RData")

## Bonferroni -----
mat_bonf <- rep(NA, length(mus))

set.seed(1)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute ruger, e-ruger, u-ruger
  rug_cr   <- apply(mat_pvals, 1, function(x) f_ruger(x, k = 1))
  rug_cr   <- sapply(rug_cr, function(x) min(x, 1))
  
  # power of a test alpha
  mat_bonf[j] <- mean(I(rug_cr <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_bonf, file = "mat_bonf2.Rdata")



