######################             Simulations           ########################
######################                Mean               ########################
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

simu_main1 <- function(data, pstar, theta0, trace = TRUE, user_lambda = NULL,
                       alpha = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                 0.6, 0.7, 0.8, 0.9, 0.95, 0.99), Sigma = NULL) {
  if(trace) print(data$id)
  y <- data$y
  n <- nrow(y)
  p <- ncol(y)

  if(is.null(Sigma)) Sigma <- diag(p)

  thetahat <- apply(y, 2, mean)

  theta_mle_oracle <- rep(0, p)
  theta_mle_oracle[1:pstar] <- thetahat[1:pstar]

  rmse_oracle <- sqrt(mean((theta0 - theta_mle_oracle)^2))
  rmse_mle <- sqrt(mean((theta0 - thetahat)^2))

  len_path <- length(user_lambda)
  if (!len_path) len_path <- 100

  scorehat <- y - matrix(rep(thetahat, each = n), nrow = n)
  # account for 1/sigma2_j in log-lik
  var.scorehat <- apply(scorehat, 2, var)
  var.scorehat.mat <- matrix(rep(var.scorehat, each = n), nrow = n)
  scorehat <- scorehat/var.scorehat.mat

  ### lambdas
  lambda.hat <- exp(seq(log(0.75), log(100), length.out = 10))

  w.hat <- matrix(NA, nrow = length(lambda.hat), ncol = p)
  tpp.hat <- fdp.hat <- nsel.hat <- rmse.hat <- aic.hat <- bic.hat <- dim.hat <-
    aic.hat.sel <- bic.hat.sel <- rmse.hat.noadj <- rep(NA, length(lambda.hat))

  for(j in 1:length(lambda.hat)) {
    if(is.character(lambda.hat[j])) {
      tpp.hat[j] <- fdp.hat[j] <- nsel.hat[j] <- rmse.hat[j] <- NA
    }
    else {
      lam <- lambda.hat[j]
      out <- cd.alg(theta = thetahat, lambda = lam, scoremat = scorehat)
      w.hat[j, ] <- out$what
      sel <- which(out$what != 0)
      pstar.hat <- length(sel)
      nsel.hat[j] <- pstar.hat
      clas <- sapply(sel, function(x) any((1:pstar) == x))
      tpp.hat[j] <- sum(clas == TRUE)/pstar
      if(length(sel) == 0) fdp.hat[j] <- 0
      else fdp.hat[j] <- sum(clas == FALSE)/nsel.hat[j]

      thetatilde <- rep(0, p)
      thetatilde[sel] <- thetahat[sel]
      rmse.hat.noadj[j] <- sqrt(mean((theta0 - thetatilde)^2))
    }
  }

  if(is.null(user_lambda))  lambda_seq <- seq(lambda.hat[1],
                            lambda.hat[length(lambda.hat)], length.out = 100)
  else lambda_seq <- user_lambda

  rmse <- nsel <- rmse_noadj <- rep(NA, length(lambda_seq))

  for(lambda in lambda_seq) {
    out <- cd.alg(theta = thetahat, lambda = lambda, scoremat = scorehat)
    sel <- which(out$what != 0)
    pstar.hat <- length(sel)
    nsel[which(lambda_seq == lambda)] <- length(sel)
    
    thetatilde <- rep(0, p)
    thetatilde[sel] <- thetahat[sel]
    rmse_noadj[which(lambda_seq == lambda)] <- sqrt(mean((theta0 - thetatilde)^2))
  }

  list(tpphat = tpp.hat, fdphat = fdp.hat, nselhat = nsel.hat, what = w.hat,
       rmsehat = rmse.hat, aichat = aic.hat, bichat = bic.hat, dimhat = dim.hat,
       lambda_seq = lambda_seq, thetahat = thetahat, nselhat = nsel.hat,
       lambdahat = lambda.hat, rmse = rmse, rmse_oracle = rmse_oracle,
       nsel = nsel, alphas = alpha, aichatsel = aic.hat.sel, bichatsel = bic.hat.sel,
       rmse_mle = rmse_mle, rmse_noadj = rmse_noadj, rmsehat_noadj = rmse.hat.noadj)
}

simu_par1 <- function(data_all, cores = 1, trace = TRUE, user_lambda = NULL,
                      alpha = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.6, 0.7, 0.8, 0.9, 0.95, 0.99), Sigma = NULL) {
  require(plyr)
  require(doMC)
  registerDoMC(ncores)

  Nsim <- length(data_all)
  pstar <- attr(data.sim, "pstar")
  theta0 <- attr(data.sim, "mu0")

  out <- llply(.data = data_all, .fun = simu_main1, pstar = pstar, trace = trace,
               .parallel = TRUE, user_lambda = user_lambda, theta0 = theta0,
               alpha = alpha, Sigma = Sigma, .inform = TRUE)
  out
}

## Setting 1
nsimu <- 2500
ncores <- 20

# sample size
n <- 250

# dimensions
d <- p <- 100
pstar <- 25
mu0 <- c(rep(5, pstar/5), rep(4, pstar/5), rep(3, pstar/5),
         rep(2, pstar/5), rep(1, pstar/5), rep(0, p - pstar))
dens <- pstar/p

# generate Gaussian sample
data.sim <- as.list(numeric(nsimu))
set.seed(321)
for(i in 1:nsimu) {
  y <- matrix(rnorm(n*p), nrow = n) + matrix(rep(mu0, each = n), nrow = n)
  data.sim[[i]] <- list(id = i, y = as.matrix(y))
}

attr(data.sim, "n") <- n
attr(data.sim, "d") <- d
attr(data.sim, "p") <- p
attr(data.sim, "mu0") <- mu0
attr(data.sim, "pstar") <- pstar
attr(data.sim, "dens") <- dens

data.sim1 <- data.sim2 <- data.sim3 <- data.sim4 <- data.sim5 <-
  as.list(numeric(500))
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
attr(data.sim1, "mu0") <- attr(data.sim2, "mu0") <- attr(data.sim3, "mu0") <-
  attr(data.sim4, "mu0") <- attr(data.sim5, "mu0") <- mu0
attr(data.sim1, "p") <- attr(data.sim2, "p") <- attr(data.sim3, "p") <-
  attr(data.sim4, "p") <- attr(data.sim5, "p") <- p
attr(data.sim1, "pstar") <- attr(data.sim2, "pstar") <- attr(data.sim3, "pstar") <-
  attr(data.sim4, "pstar") <- attr(data.sim5, "pstar") <- pstar
attr(data.sim1, "dens") <- attr(data.sim2, "dens") <- attr(data.sim3, "dens") <-
  attr(data.sim4, "dens") <- attr(data.sim5, "dens") <- dens

# grid for lambda
user_lam <- exp(seq(-14, 8.95, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam)
res <- res1
save(res, data.sim1, file = "mt_d100n250set1_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "mt_d100n250set1_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "mt_d100n250set1_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "mt_d100n250set1_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam)
res <- c(res, res5)
save(res, data.sim, file = "mt_d100n250set1.rda")

rm(data.sim, data.sim1, data.sim2, data.sim3, data.sim4, data.sim5, res, res1,
   res2, res3, res4, res5)

## Setting 2
# dimensions
d <- p <- 100
pstar <- 25
mu0 <- c(rep(5, pstar/5), rep(4, pstar/5), rep(3, pstar/5),
         rep(2, pstar/5), rep(1, pstar/5), rep(0, p - pstar))
dens <- pstar/p

### Sigma0 with rho0
rho0 <- 0.5
Sigma0 <- matrix(rho0, ncol = d, nrow = d)
diag(Sigma0) <- 1

data.sim <- as.list(numeric(nsimu))

# generate Gaussian sample
set.seed(123)
for(i in 1:nsimu) {
  y <- mvrnorm(n = n, mu = mu0, Sigma = Sigma0)
  data.sim[[i]] <- list(id = i, y = as.matrix(y))
}

attr(data.sim, "rho0") <- rho0
attr(data.sim, "n") <- n
attr(data.sim, "d") <- d
attr(data.sim, "p") <- p
attr(data.sim, "mu0") <- mu0
attr(data.sim, "pstar") <- pstar
attr(data.sim, "dens") <- dens

data.sim1 <- data.sim2 <- data.sim3 <- data.sim4 <- data.sim5 <-
  as.list(numeric(500))
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
attr(data.sim1, "mu0") <- attr(data.sim2, "mu0") <- attr(data.sim3, "mu0") <-
  attr(data.sim4, "mu0") <- attr(data.sim5, "mu0") <- mu0
attr(data.sim1, "p") <- attr(data.sim2, "p") <- attr(data.sim3, "p") <-
  attr(data.sim4, "p") <- attr(data.sim5, "p") <- p
attr(data.sim1, "pstar") <- attr(data.sim2, "pstar") <- attr(data.sim3, "pstar") <-
  attr(data.sim4, "pstar") <- attr(data.sim5, "pstar") <- pstar
attr(data.sim1, "dens") <- attr(data.sim2, "dens") <- attr(data.sim3, "dens") <-
  attr(data.sim4, "dens") <- attr(data.sim5, "dens") <- dens
attr(data.sim1, "rho0") <- attr(data.sim2, "rho0") <- attr(data.sim3, "rho0") <-
  attr(data.sim4, "rho0") <- attr(data.sim5, "rho0") <- rho0

# grid for lambda
user_lam <- exp(seq(-11, 8.95, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- res1
save(res, data.sim1, file = "mt_d100n250set2_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "mt_d100n250set2_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "mt_d100n250set2_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "mt_d100n250set2_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res5)
save(res, data.sim, file = "mt_d100n250set2.rda")