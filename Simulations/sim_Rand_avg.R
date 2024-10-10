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
mus   <- seq(0, 3, by = 0.02) # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.9                  # corr

## Randomized average -----
mat_avg <- matrix(NA, nrow = length(mus), ncol = 4)

set.seed(1234)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute f_ua, f'_ua, f*_ua
  f_ua <- matrix(NA, nrow = B, ncol = 4)
  for(i in 1:B){
    u        <- runif(1)
    f_ua[i,1] <- f_u_bavg(mat_pvals[i,], u = u)
    f_ua[i,2] <- f_u_avg(mat_pvals[i,], u = u)
    f_ua[i,3] <- f_u_avgstar(mat_pvals[i,], u = u)
    f_ua[i,4] <- f_bonf(mat_pvals[i,])
  }

  # power of a test alpha
  mat_avg[j,1] <- mean(I(f_ua[,1] <= alpha))
  mat_avg[j,2] <- mean(I(f_ua[,2] <= alpha))
  mat_avg[j,3] <- mean(I(f_ua[,3] <= alpha))
  mat_avg[j,4] <- mean(I(f_ua[,4] <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_avg, file = "mat_avg_r1.RData")

# Set hyperparameters (different rho) -----
mus   <- seq(0, 3, by = 0.02) # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.1                  # corr

## Randomized average -----
mat_avg <- matrix(NA, nrow = length(mus), ncol = 4)

set.seed(1234)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute f_ua, f'_ua, f*_ua
  f_ua <- matrix(NA, nrow = B, ncol = 4)
  for(i in 1:B){
    u        <- runif(1)
    f_ua[i,1] <- f_u_bavg(mat_pvals[i,], u = u)
    f_ua[i,2] <- f_u_avg(mat_pvals[i,], u = u)
    f_ua[i,3] <- f_u_avgstar(mat_pvals[i,], u = u)
    f_ua[i,4] <- f_bonf(mat_pvals[i,])
  }
  
  # power of a test alpha
  mat_avg[j,1] <- mean(I(f_ua[,1] <= alpha))
  mat_avg[j,2] <- mean(I(f_ua[,2] <= alpha))
  mat_avg[j,3] <- mean(I(f_ua[,3] <= alpha))
  mat_avg[j,4] <- mean(I(f_ua[,4] <= alpha))
  
  cat("iter:",j,"\n")
}

save(mat_avg, file = "mat_avg_r2.RData")

