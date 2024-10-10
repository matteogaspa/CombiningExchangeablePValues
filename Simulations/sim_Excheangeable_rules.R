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

# Set hyperparameters (different rho) -----
mus   <- seq(0, 3, by = 0.1)  # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.9                  # corr

## Exchangeable rules -----
mat_exc1 <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(123)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute exch. comb. rule
  rug_e_cr <- apply(mat_pvals, 1, function(x) f_e_ruger(x, k = K/2))
  rug_e_cr <- sapply(rug_e_cr, function(x) min(x, 1))
  
  avg_e_cr <- apply(mat_pvals, 1, function(x) f_e_bavg(x))
  avg_e_cr <- sapply(avg_e_cr, function(x) min(x, 1))
  
  harm_e_cr <- apply(mat_pvals, 1, function(x) f_e_bharm(x))
  harm_e_cr <- sapply(harm_e_cr, function(x) min(x, 1))
  
  geom_e_cr <- apply(mat_pvals, 1, function(x) f_e_bgeom(x))
  geom_e_cr <- sapply(geom_e_cr, function(x) min(x, 1))
  
  hom_e_cr  <- apply(mat_pvals, 1, function(x) f_e_hom(x, M = 50))
  hom_e_cr  <- sapply(hom_e_cr, function(x) min(x, 1))
  
  bonf_f    <- apply(mat_pvals, 1, function(x) f_bonf(x))
  bonf_f    <- sapply(bonf_f, function(x) min(x, 1))
  
  # power of a test alpha
  mat_exc1[j,1] <- mean(I(rug_e_cr <= alpha))
  mat_exc1[j,2] <- mean(I(avg_e_cr <= alpha))
  mat_exc1[j,3] <- mean(I(harm_e_cr <= alpha))
  mat_exc1[j,4] <- mean(I(geom_e_cr <= alpha))
  mat_exc1[j,5] <- mean(I(hom_e_cr <= alpha))
  mat_exc1[j,6] <- mean(I(bonf_f <= alpha))
  
  
  cat("iter:",j,"\n")
}

save(mat_exc1, file = "mat_exc1.RData")

# Set hyperparameters -----
mus   <- seq(0, 3, by = 0.1)  # vector of mu
K     <- 100                  # number of p-values
alpha <- 0.05                 # alpha level
B     <- 10^4                 # number of sims (for each mu)
rho   <- 0.1                  # corr

## Exchangeable rules -----
mat_exc2 <- matrix(NA, nrow = length(mus), ncol = 6)

set.seed(123)
for(j in 1:length(mus)){
  mu <- mus[j]
  mat_pvals <- sim_exc_pvals(n = B, K = K, rho = rho, mu = mu)
  
  # compute exch. comb. rule
  rug_e_cr <- apply(mat_pvals, 1, function(x) f_e_ruger(x, k = K/2))
  rug_e_cr <- sapply(rug_e_cr, function(x) min(x, 1))
  
  avg_e_cr <- apply(mat_pvals, 1, function(x) f_e_bavg(x))
  avg_e_cr <- sapply(avg_e_cr, function(x) min(x, 1))
  
  harm_e_cr <- apply(mat_pvals, 1, function(x) f_e_bharm(x))
  harm_e_cr <- sapply(harm_e_cr, function(x) min(x, 1))
  
  geom_e_cr <- apply(mat_pvals, 1, function(x) f_e_bgeom(x))
  geom_e_cr <- sapply(geom_e_cr, function(x) min(x, 1))
  
  hom_e_cr  <- apply(mat_pvals, 1, function(x) f_e_hom(x, M = 50))
  hom_e_cr  <- sapply(hom_e_cr, function(x) min(x, 1))
  
  bonf_f    <- apply(mat_pvals, 1, function(x) f_bonf(x))
  bonf_f    <- sapply(bonf_f, function(x) min(x, 1))
  
  # power of a test alpha
  mat_exc2[j,1] <- mean(I(rug_e_cr <= alpha))
  mat_exc2[j,2] <- mean(I(avg_e_cr <= alpha))
  mat_exc2[j,3] <- mean(I(harm_e_cr <= alpha))
  mat_exc2[j,4] <- mean(I(geom_e_cr <= alpha))
  mat_exc2[j,5] <- mean(I(hom_e_cr <= alpha))
  mat_exc2[j,6] <- mean(I(bonf_f <= alpha))
  
  
  cat("iter:",j,"\n")
}

save(mat_exc2, file = "mat_exc2.RData")

