library(psych)

# General function ------
f_calib <- function(p, M, f){
  # p: vector of p-values
  # M: number of iter
  # f: calibrator
  L <- 0
  R <- 1
  K <- length(p)
  for(m in 1:M){
    alpha <- (L+R)/2
    F_K   <- mean(f(p/alpha))
    if(F_K>=1){
      R <- alpha
    }
    else{
      L <- alpha
    }
  }
  return(R)
}

fe_calib <- function(p, M, f){
  # p: vector of p-values
  # M: number of iter
  # f: calibrator
  L <- 0
  R <- 1
  K <- length(p)
  for(m in 1:M){
    alpha <- (L+R)/2
    F_K   <- rep(NA, K)
    for(i in 1:K){
      F_K[i] <- (1/i)*sum(f(p[1:i]/alpha))
    }
    if(max(F_K)>=1){
      R <- alpha
    }
    else{
      L <- alpha
    }
  }
  return(R)
}

fu_calib <- function(p, M, f, u){
  # p: vector of p-values
  # f: calibrator
  # u: realization from a (super-)uniform r.v.
  # M: number of iter
  L <- 0
  R <- 1
  K <- length(p)
  for(m in 1:M){
    alpha <- (L+R)/2
    F_K   <- mean(f(p/alpha))
    if(F_K>=u){
      R <- alpha
    }
    else{
      L <- alpha
    }
  }
  return(R)
}



# Ruger combination rule -----
f_ruger <- function(p, k){
  # p: vector of p-values
  # k: value of k 
  K <- length(p)
  p <- p[order(p)]
  return((K/k)*p[k])
}

f_e_ruger <- function(p, k){
  # p: vector of p-values
  # k: value of k
  K <- length(p)
  comb <- rep(NA, K)
  for(i in 1:K){
    new_p   <- p[1:i]
    new_p   <- new_p[order(new_p)]
    ind     <- ceiling(i*(k/K))
    comb[i] <- (K/k)*new_p[ind]
  }
  return(min(comb))
}

f_u_ruger <- function(p, k, u){
  # p: vector of p-values
  # k: value of k
  # u: realization from a (super-)uniform r.v.
  K   <- length(p)
  p   <- p[order(p)]
  ind <- ceiling(k*u)
  return((K/k)*p[ind])
}

# "twice" the average combination rule -----
f_avg <- function(p){
  # p: vector of p-values
  return(2*mean(p))
}

f_e_avg <- function(p){
  # p: vector of p-values 
  K    <- length(p)
  comb <- rep(NA, K)
  for(i in 1:K){
    new_p   <- p[1:i]
    comb[i] <- 2*mean(new_p)
  }
  return(min(comb))
}

f_u_avg <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  return((2/(2-u))*mean(p))
}

# Improved average combination rule -----
f_bavg <- function(p){
  # p: vector of p-values
  K    <- length(p)
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(i in 1:K){
    comb[i] <- (2/max(c(2 - K/i, 0)))*mean(p[1:i])
  }
  return(min(comb))
}

f_e_bavg <- function(p){
  # p: vector of p-values
  K    <- length(p)
  mins <- rep(NA, K)
  for(l in 1:K){
    new_p <- p[1:l]
    comb  <- rep(Inf,l)
    for(m in floor(l/2):l){
      new_p   <- new_p[order(new_p)]
      comb[m] <- (2/max(c(2 - l/m, 0))) * mean(new_p[1:m])
    }
    mins[l] <- min(comb)
  }
  return(min(mins))
}

f_u_bavg <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  K    <- length(p)
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(i in 1:K){
    comb[i] <- (2/max(c(2 - (u*K)/i, 0)))*mean(p[1:i])
  }
  return(min(comb))
}

# Harmonic combination rule -----
f_harm <- function(p){
  # p: vector of p-values
  K  <- length(p)
  Tr <- log(K) + log(log(K)) + 1
  return((Tr + 1)*harmonic.mean(p))
}

f_e_harm <- function(p){
  # p: vector of p-values 
  K    <- length(p)
  Tr   <- log(K) + log(log(K)) + 1
  comb <- rep(NA, K)
  for(i in 1:K){
    new_p   <- p[1:i]
    comb[i] <- (Tr + 1)*harmonic.mean(new_p)
  }
  return(min(comb))
}

f_u_harm <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  K  <- length(p)
  Tr <- log(K) + log(log(K)) + 1
  return((u*Tr + 1)*harmonic.mean(p))
}

# Improved harmonic mean -----
f_bharm <- function(p){
  # p: vector of p-values
  K    <- length(p)
  Tr   <- log(K) + log(log(K)) + 1
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(i in 1:K){
    comb[i] <- ((K*Tr)/i + 1)*harmonic.mean(p[1:i])
  }
  return(min(comb))
}

f_e_bharm <- function(p){
  # p: vector of p-values
  K    <- length(p)
  Tr   <- log(K) + log(log(K)) + 1
  mins <- rep(NA, K)
  for(l in 1:K){
    new_p <- p[1:l]
    comb  <- rep(Inf,l)
    for(m in 1:l){
      new_p   <- new_p[order(new_p)]
      comb[m] <- ((l*Tr)/m + 1) * harmonic.mean(new_p[1:m])
    }
    mins[l] <- min(comb)
  }
  return(min(mins))
}

f_u_bharm <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  K    <- length(p)
  Tr   <- log(K) + log(log(K)) + 1
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(i in 1:K){
    comb[i] <- ((u*K*Tr)/i + 1)*harmonic.mean(p[1:i])
  }
  return(min(comb))
}

# Geometric combination rule -----
f_geom <- function(p){
  # p: vector of p-values
  return(exp(1) * geometric.mean(p))
}

f_e_geom <- function(p){
  # p: vector of p-values
  K    <- length(p)
  comb <- rep(NA, K)
  for(i in 1:K){
    new_p   <- p[1:i]
    comb[i] <- exp(1) * geometric.mean(new_p)
  }
  return(min(comb))
}

f_u_geom <- function(p, u){
  # p: vector of p-values
  return(exp(u) * geometric.mean(p))
}

# Improved geometric mean -----
f_bgeom <- function(p){
  # p: vector of p-values
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(i in 1:K){
    comb[i] <- exp(K/i) * geometric.mean(p[1:i])
  }
  return(min(comb))
}

f_e_bgeom <- function(p){
  # p: vector of p-values
  K    <- length(p)
  mins <- rep(NA, K)
  for(l in 1:K){
    new_p <- p[1:l]
    comb  <- rep(Inf,l)
    for(m in 1:l){
      new_p   <- new_p[order(new_p)]
      comb[m] <- exp(l/m) * geometric.mean(new_p[1:m])
    }
    mins[l] <- min(comb)
  }
  return(min(mins))
}

f_u_bgeom <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  K    <- length(p)
  p    <- p[order(p)]
  comb <- rep(NA, K)
  for(m in 1:K){
    comb[m] <- exp(u*(K/m)) * geometric.mean(p[1:m])
  }
  return(min(comb))
}

# Hommel function -----
grid_harm_calib <- function(p, K){
  # p: value of p
  # K: value of K
  h_K <- sum(1/(1:K))
  f   <- (K * I(h_K * p <= 1))/(ceiling(K * h_K * p))
  return(f)
}

f_e_hom <- function(p, M){
  # p: vector of p-values
  # M: number of iter
  K <- length(p)
  f <- fe_calib(p = p, M = M, f = function(x) grid_harm_calib(x, K = length(p)))
  return(f)
}

f_u_hom <- function(p, M, u){
  # p: vector of p-values
  # M: number of iter
  # u: realization from a (super-)uniform r.v.
  K <- length(p)
  f <- fu_calib(p = p, M = M, f = function(x) grid_harm_calib(x, K = length(p)), u = u)
  return(f)
}


# Randomized average introduced in Wang ----
f_u_avgstar <- function(p, u){
  # p: vector of p-values
  # u: realization from a (super-)uniform r.v.
  return(mean(p)/(2 - 2*u))
}


# Bonferroni function ----
f_bonf <- function(p){
  # p: vector of p-values
  K <- length(p)
  return(K*min(p))
}