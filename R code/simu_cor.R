######################              Simulations          ########################
# ######################        Correlation matrix       ########################
library("MASS")
library("matrixcalc")
library("mvtnorm")

#### Coordinate descent algorithm
cd.alg <- function(theta, lambda = 2, scoremat, thresh = 1e-6,
                   wstart = rep(1, length(theta)), max.it = 30) {
  p <- length(theta)
  wseq <- matrix(NA, nrow = max.it, ncol = p)
  w <- wstart
  for(i in 1:max.it) {
    for(j in 1:p) {
      uj <- scoremat[, j]
      uw_j <- apply(scoremat[, -j], 1, function(x) sum(x * w[-j]))
      numj <- sum(uj * (uj - uw_j))
      denj <- sum(uj^2)
      condit <- abs(numj) - lambda/(theta[j]^2)
      w[j] <- sign(numj) * ifelse(condit > 0, condit, 0)/denj
    }
    wseq[i, ] <- w
    if(i > 1 & (sum(wseq[i - 1, ]^2) > 0)) {
      if((sum((w - wseq[i-1, ])^2)/sum(wseq[i - 1, ]^2) < thresh)) break
    }
    if(i > 1 & (sum(wseq[i - 1, ]^2) == 0)) {
      if((sum((w - wseq[i - 1, ])^2) < thresh)) break
    }
  }
  return(list(what = w, wseq = wseq, iter = i))
}

### Bivariate marginal negative log-lik
nloglik_rho <- function(rho_st, phi_st) {
  phi_st[1]/(2*(1 - rho_st^2)) - rho_st * phi_st[2]/(1 - rho_st^2) +
    0.5 * log(1 - rho_st^2)
}

# Partial score wrt to rho_st - u_st
score_rho <- function(rho_st, phi_st) {
  - rho_st * phi_st[1]/(1 - rho_st^2)^2 + (1 + rho_st^2) *
    phi_st[2]/(1 - rho_st^2)^2 + rho_st/(1 - rho_st^2)
}

# d = 2 dimension of the sufficient statistic for rho -> ubar for i-th unit
score_bar_rho <- function(rho, phi_i, d_rho = 2) {
  p <- length(rho)
  ubar <- rep(NA, p)
  for(j in 1:p) ubar[j] <- score_rho(rho[j],
                                     phi_i[((j - 1) * d_rho + 1):(j * d_rho)])
  return(ubar)
}

### reparameterization for unconstrained optimization
beta <- function(rho) {
  0.5 * log((1 + rho)/(1 - rho))
}

rho <- function(beta) {
  (exp(2*beta) - 1)/(exp(2*beta) + 1)
}

nloglik_beta <- function(beta_st, phi_st) {
  rho_st <- rho(beta_st)
  nloglik_rho(rho_st, phi_st)
}

simu_main1 <- function(data, Sigma0, edges, true_edg, trace = TRUE,
                       user_lambda = NULL, alpha = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                                     0.6, 0.7, 0.8, 0.9, 0.95, 0.99)) {
  
  if(trace) print(data$id)
  y <- data$y
  n <- nrow(y)
  d <- ncol(y)
  p <- d * (d - 1)/2   # potential edges
  
  # standardize data
  y.sc <- scale(y)

  S <- cor(y)
  sing_flag <- is.singular.matrix(S)
  
  Sigma_mle_oracle <- matrix(0, nrow = d, ncol = d)
  Sigma_mle_oracle[Sigma0 != 0] <- S[Sigma0 != 0]
  A_mle_oracle <- Sigma_mle_oracle - Sigma0
  Delta_mle_oracle <- sqrt(sum(A_mle_oracle^2))/d
  
  A_mle <- S - Sigma0
  Delta_mle <- sqrt(sum(A_mle^2))/d
  
  len_path <- length(user_lambda)
  if (!len_path) len_path <- 100
  
  ### vector of sufficient statistics phi
  phi <- c()
  for(s in 1:(d - 1)) {
    for(t in (s + 1):d) {
      phi_st <- c(sum(y.sc[, s]^2 + y.sc[, t]^2), sum(y.sc[, s]*y.sc[, t]))
      phi <- c(phi, phi_st)
    }
  }
  
  d_rho <- 2
  ### matrix of sufficient statistics phi (n x p) for rho -> units per row
  phi_mat_rho <- matrix(NA, nrow = n, ncol = p * d_rho)
  j <- 1
  for(s in 1:(d - 1)) {
    for(t in (s + 1):d) {
      phi_st_i_rho <- apply(y.sc, 1, function(x) c(x[s]^2 + x[t]^2, x[s]*x[t]))
      phi_mat_rho[, ((j - 1)*d_rho + 1):(j*d_rho)] <- t(phi_st_i_rho)
      j <- j + 1
    }
  }
  
  ### estimate parameter rhohat from pairwise log-likelihoods
  rho0 <- 0
  
  rhohat_pw <- c()
  for(j in 1:p) {
    phi_pw <- phi[((j - 1)*d_rho + 1):(j*d_rho)]
    opt_pw <- nlminb(start = 0, objective = nloglik_beta, phi_st = phi_pw)
    rhohat_pw <- c(rhohat_pw, rho(opt_pw$par))
  }
  
  # i-th ubar by row - n x p
  # evaluated at rhohat
  scorehat <- t(apply(phi_mat_rho, 1, function(x) score_bar_rho(rho = rhohat_pw,
                                                                phi_i = x)))
  ### lambdas
  lambda.hat <- exp(seq(log(0.3), log(7), length.out = 10))
  
  w.hat <- matrix(NA, nrow = length(lambda.hat), ncol = p)
  tpp.hat <- fdp.hat <- nsel.hat <- Delta.hat <- aic.hat <- bic.hat <- sing.hat <-
    rep(NA, length(lambda.hat))
  
  for(j in 1:length(lambda.hat)) {
    if(is.character(lambda.hat[j])) {
      tpp.hat[j] <- fdp.hat[j] <- nsel.hat[j] <- Delta.hat[j] <- aic.hat[j] <-
        bic.hat[j] <- sing.hat[j] <-NA
    }
    else {
      lam <- lambda.hat[j]
      out <- cd.alg(theta = rhohat_pw, lambda = lam, scoremat = scorehat)
      w.hat[j, ] <- out$what
      sel <- which(out$what != 0)
      pstar.hat <- length(sel)
      nsel.hat[j] <- pstar.hat
      clas <- sapply(sel, function(x) any(unname(true_edg) == x))
      
      tpp.hat[j] <- sum(clas == TRUE)/length(true_edg)
      if(length(sel) == 0) fdp.hat[j] <- 0
      else fdp.hat[j] <- sum(clas == FALSE)/nsel.hat[j]
      
      rhotilde <- rep(0, p)
      rhotilde[sel] <- rhohat_pw[sel]
      
      Sigmatilde <- matrix(1, nrow = d, ncol = d)
      ind <- 1
      for(s in 1:(d - 1)) {
        for(t in (s + 1):d) {
          Sigmatilde[s, t] <- Sigmatilde[t, s] <- rhotilde[ind]
          ind <- ind + 1
        }
      }
      
      Ahat <- Sigmatilde - Sigma0
      Delta.hat[j] <- sqrt(sum(Ahat^2))/d
      sing.hat[j] <- is.singular.matrix(Sigmatilde)
    }
  }
  
  if(is.null(user_lambda))  lambda_seq <- seq(lambda.hat[1],
                                              lambda.hat[length(lambda.hat)], length.out = 100)
  else lambda_seq <- user_lambda
  
  Delta <- nsel <- Delta2 <- rep(NA, length(lambda_seq))
  
  for(lambda in lambda_seq) {
    out <- cd.alg(theta = rhohat_pw, lambda = lambda, scoremat = scorehat)
    sel <- which(out$what != 0)
    pstar.hat <- length(sel)
    nsel[which(lambda_seq == lambda)] <- length(sel)
    
    rhotilde <- rep(0, p)
    rhotilde[sel] <- rhohat_pw[sel]
    
    Sigmatilde <- matrix(1, nrow = d, ncol = d)
    ind <- 1
    for(s in 1:(d - 1)) {
      for(t in (s + 1):d) {
        Sigmatilde[s, t] <- Sigmatilde[t, s] <- rhotilde[ind]
        ind <- ind + 1
      }
    }
    
    A <- Sigmatilde - Sigma0
    Delta[which(lambda_seq == lambda)] <- sqrt(sum(A^2))/d
    diff <- rhotilde - as.vector(Sigma0[lower.tri(Sigma0)])
    Delta2[which(lambda_seq == lambda)] <- sqrt(sum(diff^2))
  }
  
  list(tpphat = tpp.hat, fdphat = fdp.hat, nselhat = nsel.hat, what = w.hat,
       Deltahat = Delta.hat, aichat = aic.hat, bichat = bic.hat, singS = sing_flag,
       lambda_seq = lambda_seq, S = S, singhat = sing.hat, Delta2 = Delta2,
       lambdahat = lambda.hat, Delta = Delta, nsel = nsel, alphas = alpha,
       Delta_oracle = Delta_mle_oracle, Delta_mle = Delta_mle)
}

simu_par1 <- function(data_all, cores = 1, trace = TRUE, user_lambda = NULL,
                      alpha = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 0.95, 0.99)) {
  require(plyr)
  require(doMC)
  registerDoMC(ncores)
  
  Nsim <- length(data_all)
  Sigma0 <- attr(data_all, "Sigma0")
  true_edges <- attr(data.sim, "true_edges")
  edges <- attr(data.sim, "edges")
  
  out <- llply(.data = data_all, .fun = simu_main1, Sigma0 = Sigma0,
               true_edg = true_edges, edges = edges, trace = trace,
               .parallel = TRUE, user_lambda = user_lambda, .inform = TRUE,
               alpha = alpha)
  out
}

## Setting 1
nsimu <- 2500
ncores <- 20

## sample size
nobs <- 250

### Sigma0
d <- 15
p <- d * (d - 1)/2

rho0 <- 0.5

###########################################
Sigma0 <- matrix(rep(0, d * d), ncol = d)
Sigma0[1:5, 1:5] <- rho0
diag(Sigma0) <- 1
###########################################

# edges configuration
edges <- rep(NA, p)
i <- 1
for(s in 1:(d - 1)) {
  for(t in (s + 1):d) {
    edges[i] <- Sigma0[s, t]
    names(edges)[i] <- paste0("(", s, ",", t, ")")
    i <- i + 1
  }
}

true_edges <- which(edges != 0)    # nonzero edges
n_edg <- length(true_edges)

dens_mat <- n_edg/p
dens_mat

data.sim <- as.list(numeric(nsimu))

# generate Gaussian sample
set.seed(123)
for(i in 1:nsimu) {
  y <- rmvnorm(n = nobs, mean = rep(0, d), sigma = Sigma0)
  data.sim[[i]] <- list(id = i, y = as.matrix(y))
}

attr(data.sim, "n") <- nobs
attr(data.sim, "d") <- d
attr(data.sim, "rho0") <- rho0
attr(data.sim, "Sigma0") <- Sigma0
attr(data.sim, "true_edges") <- true_edges
attr(data.sim, "n_edges") <- n_edg
attr(data.sim, "edges") <- edges
attr(data.sim, "dens_mat") <- dens_mat
attr(data.sim, "n_cliques") <- ncliq

data.sim1 <- data.sim2 <- data.sim3 <- data.sim4 <- data.sim5 <-
  as.list(numeric(ncores))
for(i in 1:500) {
  data.sim1[[i]] <- data.sim[[i]]
  data.sim2[[i]] <- data.sim[[500 + i]]
  data.sim3[[i]] <- data.sim[[1000 + i]]
  data.sim4[[i]] <- data.sim[[1500 + i]]
  data.sim5[[i]] <- data.sim[[2000 + i]]
}

attr(data.sim1, "n") <- attr(data.sim2, "n") <- attr(data.sim3, "n") <-
  attr(data.sim4, "n") <- attr(data.sim5, "n") <- nobs
attr(data.sim1, "d") <- attr(data.sim2, "d") <- attr(data.sim3, "d") <-
  attr(data.sim4, "d") <- attr(data.sim5, "d") <- d
attr(data.sim1, "rho0") <- attr(data.sim2, "rho0") <- attr(data.sim3, "rho0") <-
  attr(data.sim4, "rho0") <- attr(data.sim5, "rho0") <- rho0
attr(data.sim1, "Sigma0") <- attr(data.sim2, "Sigma0") <-
  attr(data.sim3, "Sigma0") <-
  attr(data.sim4, "Sigma0") <- attr(data.sim5, "Sigma0") <- Sigma0
attr(data.sim1, "true_edges") <- attr(data.sim2, "true_edges") <-
  attr(data.sim3, "true_edges") <-
  attr(data.sim4, "true_edges") <- attr(data.sim5, "true_edges") <- true_edges
attr(data.sim1, "n_edges") <- attr(data.sim2, "n_edges") <-
  attr(data.sim3, "n_edges") <-
  attr(data.sim4, "n_edges") <- attr(data.sim5, "n_edges") <- n_edg
attr(data.sim1, "edges") <- attr(data.sim2, "edges") <-
  attr(data.sim3, "edges") <-
  attr(data.sim4, "edges") <- attr(data.sim5, "edges") <- edges
attr(data.sim1, "dens_mat") <- attr(data.sim2, "dens_mat") <-
  attr(data.sim3, "dens_mat") <-
  attr(data.sim4, "dens_mat") <- attr(data.sim5, "dens_mat") <- dens_mat
attr(data.sim1, "n_cliques") <- attr(data.sim2, "n_cliques") <-
  attr(data.sim3, "n_cliques") <-
  attr(data.sim4, "n_cliques") <- attr(data.sim5, "n_cliques") <- ncliq

user_lam <- exp(seq(-15, 6, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam)
res <- res1
save(res, data.sim1, file = "cs_d15n250set1_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "cs_d15n250set1_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "cs_d15n250set1_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "cs_d15n250set1_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam)
res <- c(res, res5)
save(res, data.sim, file = "cs_d15n250set1.rda")

## Setting 2 - Toeplitz (delta = 0.1)
# sample size
nobs <- 250

### Sigma0
d <- 15
p <- d * (d - 1)/2    # total number of edges

#########################################
delta <- 0.1
Sigma0 <- Sigmatot <- matrix(rep(0, d * d), ncol = d)
for(s in 1:d) {
  for(t in 1:d) Sigmatot[s, t] <- exp(-delta * abs(t - s))
}
Sigma0[1:5, 1:5] <- Sigmatot[1:5, 1:5]
diag(Sigma0) <- 1
#########################################

# edges configuration
edges <- rep(NA, p)
i <- 1
for(s in 1:(d - 1)) {
  for(t in (s + 1):d) {
    edges[i] <- Sigma0[s, t]
    names(edges)[i] <- paste0("(", s, ",", t, ")")
    i <- i + 1
  }
}

true_edges <- which(edges != 0)    # nonzero edges
n_edg <- length(true_edges)

dens_mat <- n_edg/p
dens_mat

data.sim <- as.list(numeric(nsimu))

# generate Gaussian sample
set.seed(123)
for(i in 1:nsimu) {
  y <- rmvnorm(n = nobs, mean = rep(0, d), sigma = Sigma0)
  data.sim[[i]] <- list(id = i, y = as.matrix(y))
}

attr(data.sim, "n") <- nobs
attr(data.sim, "d") <- d
attr(data.sim, "Sigma0") <- Sigma0
attr(data.sim, "true_edges") <- true_edges
attr(data.sim, "n_edges") <- n_edg
attr(data.sim, "edges") <- edges
attr(data.sim, "dens_mat") <- dens_mat

data.sim1 <- data.sim2 <- data.sim3 <- data.sim4 <- data.sim5 <-
  as.list(numeric(ncores))
for(i in 1:500) {
  data.sim1[[i]] <- data.sim[[i]]
  data.sim2[[i]] <- data.sim[[500 + i]]
  data.sim3[[i]] <- data.sim[[1000 + i]]
  data.sim4[[i]] <- data.sim[[1500 + i]]
  data.sim5[[i]] <- data.sim[[2000 + i]]
}

attr(data.sim1, "n") <- attr(data.sim2, "n") <- attr(data.sim3, "n") <-
  attr(data.sim4, "n") <- attr(data.sim5, "n") <- nobs
attr(data.sim1, "d") <- attr(data.sim2, "d") <- attr(data.sim3, "d") <-
  attr(data.sim4, "d") <- attr(data.sim5, "d") <- d
attr(data.sim1, "Sigma0") <- attr(data.sim2, "Sigma0") <-
  attr(data.sim3, "Sigma0") <-
  attr(data.sim4, "Sigma0") <- attr(data.sim5, "Sigma0") <- Sigma0
attr(data.sim1, "true_edges") <- attr(data.sim2, "true_edges") <-
  attr(data.sim3, "true_edges") <-
  attr(data.sim4, "true_edges") <- attr(data.sim5, "true_edges") <- true_edges
attr(data.sim1, "n_edges") <- attr(data.sim2, "n_edges") <-
  attr(data.sim3, "n_edges") <-
  attr(data.sim4, "n_edges") <- attr(data.sim5, "n_edges") <- n_edg
attr(data.sim1, "edges") <- attr(data.sim2, "edges") <-
  attr(data.sim3, "edges") <-
  attr(data.sim4, "edges") <- attr(data.sim5, "edges") <- edges
attr(data.sim1, "dens_mat") <- attr(data.sim2, "dens_mat") <-
  attr(data.sim3, "dens_mat") <-
  attr(data.sim4, "dens_mat") <- attr(data.sim5, "dens_mat") <- dens_mat

user_lam <- exp(seq(-15, 6, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam)
res <- res1
save(res, data.sim1, file = "cs_d15n250set2_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "cs_d15n250set2_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "cs_d15n250set2_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "cs_d15n250set2_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam)
res <- c(res, res5)
save(res, data.sim, file = "cs_d15n250set2.rda")
