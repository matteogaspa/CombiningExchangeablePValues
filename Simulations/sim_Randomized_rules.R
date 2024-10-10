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
mus   <- seq(0, 3, by = 0.05) # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.9                  # corr

mat_ran <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(1234)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute randomized rules
  f_u <- matrix(NA, nrow = B, ncol = 6)
  for(i in 1:B){
    u        <- runif(1)
    f_u[i,1] <- f_u_ruger(mat_pvals[i,], k = K/2, u = u)
    f_u[i,2] <- f_u_bavg(mat_pvals[i,], u = u)
    f_u[i,3] <- f_u_bgeom(mat_pvals[i,], u = u)
    f_u[i,4] <- f_u_bharm(mat_pvals[i,], u = u)
    f_u[i,5] <- f_u_hom(mat_pvals[i,], M = 100, u = u)
    f_u[i,6] <- f_bonf(mat_pvals[i,])
  }
  
  # power of a test alpha
  mat_ran[j,1] <- mean(I(f_u[,1] <= alpha))
  mat_ran[j,2] <- mean(I(f_u[,2] <= alpha))
  mat_ran[j,3] <- mean(I(f_u[,3] <= alpha))
  mat_ran[j,4] <- mean(I(f_u[,4] <= alpha))
  mat_ran[j,5] <- mean(I(f_u[,5] <= alpha))
  mat_ran[j,6] <- mean(I(f_u[,6] <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_ran, file = "mat_ran1.RData")

# Set hyperparameters (different rho) -----
mus   <- seq(0, 3, by = 0.05) # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.1                  # corr

mat_ran <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(1234)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute f_ua, f'_ua, f*_ua
  f_u <- matrix(NA, nrow = B, ncol = 6)
  for(i in 1:B){
    u        <- runif(1)
    f_u[i,1] <- f_u_ruger(mat_pvals[i,], k = K/2, u = u)
    f_u[i,2] <- f_u_bavg(mat_pvals[i,], u = u)
    f_u[i,3] <- f_u_bgeom(mat_pvals[i,], u = u)
    f_u[i,4] <- f_u_bharm(mat_pvals[i,], u = u)
    f_u[i,5] <- f_u_hom(mat_pvals[i,], M = 100, u = u)
    f_u[i,6] <- f_bonf(mat_pvals[i,])
  }
  
  # power of a test alpha
  mat_ran[j,1] <- mean(I(f_u[,1] <= alpha))
  mat_ran[j,2] <- mean(I(f_u[,2] <= alpha))
  mat_ran[j,3] <- mean(I(f_u[,3] <= alpha))
  mat_ran[j,4] <- mean(I(f_u[,4] <= alpha))
  mat_ran[j,5] <- mean(I(f_u[,5] <= alpha))
  mat_ran[j,6] <- mean(I(f_u[,6] <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_ran, file = "mat_ran2.RData")
