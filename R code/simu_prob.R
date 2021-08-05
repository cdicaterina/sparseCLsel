######################              Simulations          ########################
######################        Probit with Covariate      ########################
library("MASS")
library("matrixcalc")
library("mvtnorm")

nll.prob <- function(theta, data.y, x.cov) {
  mu <- theta[1] + theta[2]*x.cov
  -sum(data.y * log(pnorm(mu)) + (1 - data.y) * log(1 - pnorm(mu)))
}

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
  x <- data$x
  n <- nrow(y)
  p <- ncol(y)

  if(is.null(Sigma)) Sigma <- diag(p)

  betahat <- apply(y, 2, function(z) nlminb(start = c(1, 1), objective = nll.prob,
                                            data.y = z, x.cov = x)$par)
  thetahat <- betahat[2, ]

  beta0mat <- matrix(rep(betahat[1, ], each = n), nrow = n)
  beta1mat <- matrix(rep(thetahat, each = n), nrow = n)
  xmat <- matrix(rep(x, times = p), nrow = n)
  mumat <- beta0mat + beta1mat * xmat

  theta_mle_oracle <- rep(0, p)
  theta_mle_oracle[1:pstar] <- thetahat[1:pstar]

  rmse_oracle <- sqrt(mean((theta0 - theta_mle_oracle)^2))
  rmse_mle <- sqrt(mean((theta0 - thetahat)^2))

  len_path <- length(user_lambda)
  if (!len_path) len_path <- 100

  scorehat <- dnorm(mumat) * (y - pnorm(mumat)) * xmat/
    (pnorm(mumat) * (1 - pnorm(mumat)))

  ### lambdas
  lambda.hat <- exp(seq(log(0.2), log(40), length.out = 10))

  w.hat <- matrix(NA, nrow = length(lambda.hat), ncol = p)
  tpp.hat <- fdp.hat <- nsel.hat <- rmse.hat <- aic.hat <- bic.hat <- dim.hat <-
    aic.hat.sel <- bic.hat.sel <- rmse.hat.noadj <- rep(NA, length(lambda.hat))

  for(j in 1:length(lambda.hat)) {
    if(is.character(lambda.hat[j])) {
      tpp.hat[j] <- fdp.hat[j] <- nsel.hat[j] <- rmse.hat[j] <- rmse.hat.noadj[j] <- NA
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
       rmsehat = rmse.hat, dimhat = dim.hat,
       lambda_seq = lambda_seq, thetahat = thetahat, lambdahat = lambda.hat,
       rmse = rmse, rmse_oracle = rmse_oracle, nsel = nsel, alphas = alpha,
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
  theta0 <- attr(data.sim, "beta10")

  out <- llply(.data = data_all, .fun = simu_main1, pstar = pstar, trace = trace,
               .parallel = TRUE, user_lambda = user_lambda, theta0 = theta0,
               alpha = alpha, Sigma = Sigma, .inform = TRUE)
  out
}

# Setting 1
nsimu <- 2500
ncores <- 20

# sample size
n <- 250

# dimensions
d <- p <- 100
pstar <- 25
# true parameters
beta00 <- rep(0.1, times = d)
beta10 <- c(rep(1.5, pstar/5), rep(1.25, pstar/5), rep(1, pstar/5),
            rep(0.75, pstar/5), rep(0.5, pstar/5), rep(0, p - pstar))
dens <- pstar/p
# covariate
set.seed(123)
x <- rnorm(n, mean = 0)

# matrix form
beta00mat <- matrix(rep(beta00, each = n), nrow = n)
beta10mat <- matrix(rep(beta10, each = n), nrow = n)
xmat <- matrix(rep(x, times = d), nrow = n)
mumat <- beta00mat + beta10mat*xmat

# generate Gaussian sample
data.sim <- as.list(numeric(nsimu))
set.seed(321)
for(i in 1:nsimu) {
  ystar <- matrix(rnorm(n*p), nrow = n) + mumat
  y <- ifelse(ystar > 0, 1, 0)
  data.sim[[i]] <- list(id = i, y = as.matrix(y), x = x)
}

attr(data.sim, "n") <- n
attr(data.sim, "d") <- d
attr(data.sim, "p") <- p
attr(data.sim, "beta00") <- beta00
attr(data.sim, "beta10") <- beta10
attr(data.sim, "pstar") <- pstar
attr(data.sim, "dens") <- dens
attr(data.sim, "Sigma") <- diag(d)

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
attr(data.sim1, "beta00") <- attr(data.sim2, "beta00") <- attr(data.sim3, "beta00") <-
  attr(data.sim4, "beta00") <- attr(data.sim5, "beta00") <- beta00
attr(data.sim1, "beta10") <- attr(data.sim2, "beta10") <- attr(data.sim3, "beta10") <-
  attr(data.sim4, "beta10") <- attr(data.sim5, "beta10") <- beta10
attr(data.sim1, "p") <- attr(data.sim2, "p") <- attr(data.sim3, "p") <-
  attr(data.sim4, "p") <- attr(data.sim5, "p") <- p
attr(data.sim1, "pstar") <- attr(data.sim2, "pstar") <- attr(data.sim3, "pstar") <-
  attr(data.sim4, "pstar") <- attr(data.sim5, "pstar") <- pstar
attr(data.sim1, "dens") <- attr(data.sim2, "dens") <- attr(data.sim3, "dens") <-
  attr(data.sim4, "dens") <- attr(data.sim5, "dens") <- dens
attr(data.sim1, "Sigma") <- attr(data.sim2, "Sigma") <- attr(data.sim3, "Sigma") <-
  attr(data.sim4, "Sigma") <- attr(data.sim5, "Sigma") <- diag(d)

# grid for lambda
user_lam <- exp(seq(-8.5, 5, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam)
res <- res1
save(res, data.sim1, file = "prob_d100n250set1_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "prob_d100n250set1_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "prob_d100n250set1_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "prob_d100n250set1_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam)
res <- c(res, res5)
save(res, data.sim, file = "prob_d100n250set1.rda")

## Setting 2
# dimensions
d <- p <- 100
pstar <- 25

# true parameters
beta00 <- rep(0.1, times = d)
beta10 <- c(rep(1.5, pstar/5), rep(1.25, pstar/5), rep(1, pstar/5),
            rep(0.75, pstar/5), rep(0.5, pstar/5), rep(0, p - pstar))

# covariate
set.seed(123)
x <- rnorm(n, mean = 0)

# matrix form
beta00mat <- matrix(rep(beta00, each = n), nrow = n)
beta10mat <- matrix(rep(beta10, each = n), nrow = n)
xmat <- matrix(rep(x, times = d), nrow = n)
mumat <- beta00mat + beta10mat*xmat

### Sigma0 with rho0
rho0 <- 0.5
Sigma0 <- matrix(rho0, ncol = d, nrow = d)
diag(Sigma0) <- 1

data.sim <- as.list(numeric(nsimu))

# generate Gaussian sample
set.seed(123)
for(i in 1:nsimu) {
  ystar <- mvrnorm(n = n, mu = rep(0, d), Sigma = Sigma0) + mumat
  y <- ifelse(ystar > 0, 1, 0)
  data.sim[[i]] <- list(id = i, y = as.matrix(y), x = x)
}

attr(data.sim, "rho0") <- rho0
attr(data.sim, "n") <- n
attr(data.sim, "d") <- d
attr(data.sim, "p") <- p
attr(data.sim, "beta00") <- beta00
attr(data.sim, "beta10") <- beta10
attr(data.sim, "pstar") <- pstar
attr(data.sim, "dens") <- dens
attr(data.sim, "Sigma") <- Sigma0

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
attr(data.sim1, "beta00") <- attr(data.sim2, "beta00") <- attr(data.sim3, "beta00") <-
  attr(data.sim4, "beta00") <- attr(data.sim5, "beta00") <- beta00
attr(data.sim1, "beta10") <- attr(data.sim2, "beta10") <- attr(data.sim3, "beta10") <-
  attr(data.sim4, "beta10") <- attr(data.sim5, "beta10") <- beta10
attr(data.sim1, "p") <- attr(data.sim2, "p") <- attr(data.sim3, "p") <-
  attr(data.sim4, "p") <- attr(data.sim5, "p") <- p
attr(data.sim1, "pstar") <- attr(data.sim2, "pstar") <- attr(data.sim3, "pstar") <-
  attr(data.sim4, "pstar") <- attr(data.sim5, "pstar") <- pstar
attr(data.sim1, "dens") <- attr(data.sim2, "dens") <- attr(data.sim3, "dens") <-
  attr(data.sim4, "dens") <- attr(data.sim5, "dens") <- dens
attr(data.sim1, "Sigma") <- attr(data.sim2, "Sigma") <- attr(data.sim3, "Sigma") <-
  attr(data.sim4, "Sigma") <- attr(data.sim5, "Sigma") <- Sigma0

# grid for lambda
user_lam <- exp(seq(-8.5, 5, length.out = 100))
res1 <- simu_par1(data.sim1, cores = ncores, user_lambda = user_lam,
                  Sigma = Sigma0)
res <- res1
save(res, data.sim1, file = "prob_d100n250set2_1.rda")

res2 <- simu_par1(data.sim2, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res2)
data.sim12 <- c(data.sim1, data.sim2)
save(res, data.sim12, file = "prob_d100n250set2_2.rda")

res3 <- simu_par1(data.sim3, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res3)
data.sim123 <- c(data.sim1, data.sim2, data.sim3)
save(res, data.sim123, file = "prob_d100n250set2_3.rda")

res4 <- simu_par1(data.sim4, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res4)
data.sim1234 <- c(data.sim1, data.sim2, data.sim3, data.sim4)
save(res, data.sim1234, file = "prob_d100n250set2_4.rda")

res5 <- simu_par1(data.sim5, cores = ncores, user_lambda = user_lam, Sigma = Sigma0)
res <- c(res, res5)
save(res, data.sim, file = "prob_d100n250set2.rda")

