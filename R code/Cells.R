#########################################################################################
#############################          Cells data       #################################
#########################################################################################
library("readxl")

library("FinCovRegularization")
library("spcov")

setwd("~/Dropbox/DiCaterinaFerrari2020/Supplementary/sparseCLsel/cellsData")

cell1 <- read_xls("1. cd3cd28.xls")
cell2 <- read_xls("2. cd3cd28icam2.xls")
cell3 <- read_xls("3. cd3cd28+aktinhib.xls")
cell4 <- read_xls("4. cd3cd28+g0076.xls")
cell5 <- read_xls("5. cd3cd28+psitect.xls")
cell6 <- read_xls("6. cd3cd28+u0126.xls")
cell7 <- read_xls("7. cd3cd28+ly.xls")
cell8 <- read_xls("8. pma.xls")
cell9 <- read_xls("9. b2camp.xls")

names(cell8) <- names(cell7)

cell <- rbind(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9)
rm(cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9)

names(cell) <- c("Raf", "Mek", "Plcg", "PIP2", "PIP3", "Erk", "Akt",
                 "PKA", "PKC", "P38", "Jnk")

# dimensions
d <- ncol(cell)
p <- d * (d - 1)/2    # total number of edges

# edges configuration
rho_est <- sigma_est <- rhomm_est <- rep(0, p)
i <- 1
for(s in 1:(d - 1)) {
  for(t in (s + 1):d) {
    names(rho_est)[i] <- names(sigma_est)[i] <- paste0("(", dimnames(cell)[[2]][s],
                                              ",", dimnames(cell)[[2]][t], ")")
    names(rhomm_est)[i] <- names(rho_est)[i]
    i <- i + 1
  }
}

cell_mat <- as.matrix(cell)
y <- cell_mat
n <- nrow(y)

# sample correlation matrix
cormat <- cor(y)
rhohat <- as.vector(cormat[lower.tri(cormat)])
named_rhohat <- rhohat
index <- 1:p
names(index) <- names(named_rhohat) <- names(rho_est)

# order the ML estimates of rhos
orderedMLE <- sort(abs(named_rhohat), decreasing = TRUE)
mle.ind <- index[names(orderedMLE)]

k <- 6
mle.ind[1:k]
orderedMLE[1:k]

### Find SEs for ML estimates of correlation coefficients (Nel, 1985):
SEs <- sqrt(1/n * (1 + named_rhohat^4 - named_rhohat^2))
# z-stats
z.scores <- named_rhohat/SEs
sort(abs(z.scores), decreasing = TRUE)[1:6]

## Standardize data
ysd <- scale(y)

#### Coordinate descent algorithm 
cd.alg <- function(theta, lambda = 2, scoremat, thresh = 1e-6, 
                   wstart = rep(1, length(theta)), max.it = 30) {
  p <- length(theta)
  wseq <- matrix(NA, nrow = max.it, ncol = p)
  w <- wstart
  for(i in 1:max.it) {
    if(i %% 10 == 0) cat("iteration", i, '\n')
    for(j in 1:p) {
      uj <- scoremat[, j]
      uw_j <- apply(scoremat[, -j], 1, function(x) sum(x * w[-j]))
      numj <- sum(uj * (uj - uw_j))
      denj <- sum(uj^2)
      condit <- abs(numj) - lambda/(theta[j]^2)
      w[j] <- sign(numj) * ifelse(condit > 0, condit, 0)/denj
    }
    wseq[i, ] <- w
    if((i > 1) & (sum((w - wseq[i-1, ])^2)/sum(wseq[i - 1,]^2) < thresh)) break
  }
  return(list(what = w, wseq = wseq, iter = i))
}

### Pairwise likelihood approach
# Partial score wrt to rho_st - u_st
score_rho <- function(rho_st, phi_st) {
  - rho_st * phi_st[1]/(1 - rho_st^2)^2 + (1 + rho_st^2) * 
    phi_st[2]/(1 - rho_st^2)^2 + rho_st/(1 - rho_st^2)
}

# d = 2 is the dimension of the sufficient statistic for rho -> ubar for i-th unit
score_bar_rho <- function(rho, phi_i, d_rho = 2) {
  m <- length(rho)
  ubar <- rep(NA, m)
  for(j in 1:m) ubar[j] <- score_rho(rho[j], 
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

### Bivariate marginal negative log-lik
nloglik_rho <- function(rho_st, phi_st) {
  phi_st[1]/(2*(1 - rho_st^2)) - rho_st * phi_st[2]/(1 - rho_st^2) +
    0.5 * log(1 - rho_st^2)
}

nloglik_beta <- function(beta_st, phi_st) {
  rho_st <- rho(beta_st)
  nloglik_rho(rho_st, phi_st)
}

### vector of sufficient statistics phi
phi <- c()
for(s in 1:(d - 1)) {
  for(t in (s + 1):d) {
    phi_st <- c(sum(ysd[, s]^2 + ysd[, t]^2), sum(ysd[, s]*ysd[, t]))
    phi <- c(phi, phi_st)
  }
}

d_rho <- 2
### matrix of sufficient statistics phi (n x p) for rho -> units per row
phi_mat_rho <- matrix(NA, nrow = n, ncol = p * d_rho)
j <- 1
for(s in 1:(d - 1)) {
  for(t in (s + 1):d) {
    phi_st_i_rho <- apply(ysd, 1, function(x) c(x[s]^2 + x[t]^2, x[s]*x[t]))
    phi_mat_rho[, ((j - 1)*d_rho + 1):(j*d_rho)] <- t(phi_st_i_rho)
    j <- j + 1
  }
}

### estimate parameter rhohat from pairwise log-likelihoods
rhohat_pw <- c()
for(j in 1:p) {
  phi_pw <- phi[((j - 1)*d_rho + 1):(j*d_rho)]
  opt_pw <- nlminb(start = 0, objective = nloglik_beta, phi_st = phi_pw)
  rhohat_pw <- c(rhohat_pw, rho(opt_pw$par))
}

# i-th ubar by row - n x p
# evaluated at rhohat
Ubar <- t(apply(phi_mat_rho, 1, function(x) score_bar_rho(rho = rhohat_pw, 
                                                          phi_i = x)))

lambda.val <- 140
# lambda.val <- 650
# lambda.val <- 7000

## Fit
out.cell <- cd.alg(theta = rhohat_pw, lambda = lambda.val, scoremat = Ubar, 
                   max.it = 30)
w <- out.cell$what

# selection
nsel <- sum(w != 0)
nsel

rho_est <- rep(0, p)
names(rho_est) <- names(named_rhohat)
rho_est[w != 0] <- rhohat_pw[w != 0]

sel <- which(rho_est != 0)
sel

## lambda = 140
# [1] 25  
# (Raf,Mek)  (Mek,Plcg)  (Mek,PIP2)   (Mek,Akt)   (Mek,P38) (Plcg,PIP2) 
# 1          11          12          15          18          20 
# (Plcg,Akt)  (Plcg,PKA)  (Plcg,PKC)  (Plcg,P38)  (Plcg,Jnk) (PIP2,PIP3) 
# 23          24          25          26          27          28 
# (PIP2,Akt)  (PIP2,PKC)  (PIP2,P38)  (PIP2,Jnk)   (Erk,Akt)   (Erk,PKA) 
# 30          32          33          34          41          42 
# (Erk,Jnk)   (Akt,PKC)   (Akt,P38)   (Akt,Jnk)   (PKC,P38)   (PKC,Jnk) 
# 45          47          48          49          53          54 
# (P38,Jnk) 
# 55 

## lambda = 650
# [1] 12  
# (Raf,Mek) (Plcg,PIP2)  (Plcg,Akt)  (Plcg,P38)  (Plcg,Jnk)  (PIP2,Akt) 
# 1          20          23          26          27          30 
# (PIP2,P38)   (Erk,Akt)   (Akt,Jnk)   (PKC,P38)   (PKC,Jnk)   (P38,Jnk) 
# 33          41          49          53          54          55

## lambda = 7000
# [1] 6   
# (Plcg,PIP2)  (Plcg,Akt)   (Erk,Akt)   (PKC,P38)   (PKC,Jnk)   (P38,Jnk) 
# 20          23          41          53          54          55

#############################       Competing methods       ###########################
len_path <- 100
# covariance matrix
S <- cov(y)

## Soft thresholding 
cv_soft <- threshold.cv(y, method = "soft", thresh.len = len_path, n.cv = 30,
                        norm = "F", seed = 321)
thresh_soft <- cv_soft$threshold.grid

threshat_soft <- cv_soft$parameter.opt
nsel_soft <- rep(NA, length(thresh_soft))
sel_soft <- list(as.numeric(len_path))
for(j in 1:length(thresh_soft)) {
  Sigma_soft <- soft.thresholding(S, threshold = thresh_soft[j])
  edges_soft <- as.vector(Sigma_soft[lower.tri(Sigma_soft)])

  sel_soft[[j]] <- which(edges_soft != 0)
  nsel_soft[j] <- length(which(edges_soft != 0))
}

# choose pstarhat for comparison
pstarhat <- 6
ind.p <- which(nsel_soft == pstarhat)[1]
indexes <- sel_soft[[ind.p]]
rhohat_soft <- rep(0, p)
index <- 1:p
names(index) <- names(rhohat_soft) <- names(rho_est)
rhohat_soft[indexes] <- rhohat[indexes]

sel.thr <- index[indexes]
sel.thr
# (Raf,Mek)   (Mek,P38) (Plcg,PIP2)  (PIP2,P38)   (PKA,P38)   (P38,Jnk) 
# 1          18          20          33          51          55 

### Bien & Tibshirani (2011)
P <- matrix(1, d, d)                      
diag(P) <- 0

lambda_bt <- seq(0.0001, 0.1, length.out = 50)
sel_bt <- list()
nsel_bt <- rep(NA, length(lambda_bt))

if(!is.null(lambda_bt))  {
  for(j in 1:length(lambda_bt)) {
    Sigma_bt <- spcov(Sigma = S, S = S, lambda = lambda_bt[j] * P,
                      step.size = 100)$Sigma
    edges_bt <- as.vector(Sigma_bt[lower.tri(Sigma_bt)])
    sel_bt[[j]] <- which(edges_bt != 0)
    nsel_bt[j] <- length(which(edges_bt != 0))
  }
}

# choose pstarhat for comparison
pstarhat <- 6
ind.p <- which(nsel_bt == pstarhat)[1]
ind.p

sel.bt <- names(rho_est)[sel_bt[[ind.p]]]
sel.bt
# [1] "(Raf,Mek)"   "(Mek,P38)"   "(Plcg,PIP2)" "(PIP2,P38)"  "(PKA,P38)"   "(P38,Jnk)"  


################################################################################
##################                 EXPERIMENT                ###################
################################################################################
### Data split in random subsets with (approximately) nstar obs
nstar <- 250
ngroups <- ceiling(nrow(cell_mat)/nstar)

ss <- matrix(NA, ncol = ngroups, nrow = nstar)
set.seed(123)
ind <- sample(1:nrow(cell_mat), nstar, replace = F)
ss[, 1] <- ind
for(l in 2:(ngroups - 1)) {
  ss[, l] <- sample((1:nrow(cell_mat))[-ind], nstar, replace = F)
  ind <- c(ind, ss[, l])
}
ss[1:(length((1:nrow(cell_mat))[-ind])), ngroups] <- (1:nrow(cell_mat))[-ind]

rmse <- rep(NA, ngroups)
nsel <- rep(NA, ngroups)

for(i in 1:ngroups) {
  sub.i <- ss[!is.na(ss[, i]), i]
  y <- cell_mat[sub.i, ]
  n <- nrow(y)
  
  # sample correlation matrix
  cormat <- cor(y)
  rhohat <- as.vector(cormat[lower.tri(cormat)])
  named_rhohat <- rhohat
  index <- 1:p
  names(index) <- names(named_rhohat) <- names(rho_est)
  
  ## Standardize data
  ysd <- scale(y)
  
  ### vector of sufficient statistics phi
  phi <- c()
  for(s in 1:(d - 1)) {
    for(t in (s + 1):d) {
      phi_st <- c(sum(ysd[, s]^2 + ysd[, t]^2), sum(ysd[, s]*ysd[, t]))
      phi <- c(phi, phi_st)
    }
  }
  
  d_rho <- 2
  ### matrix of sufficient statistics phi (n x p) for rho -> units per row
  phi_mat_rho <- matrix(NA, nrow = n, ncol = p * d_rho)
  j <- 1
  for(s in 1:(d - 1)) {
    for(t in (s + 1):d) {
      phi_st_i_rho <- apply(ysd, 1, function(x) c(x[s]^2 + x[t]^2, x[s]*x[t]))
      phi_mat_rho[, ((j - 1)*d_rho + 1):(j*d_rho)] <- t(phi_st_i_rho)
      j <- j + 1
    }
  }
  
  ### estimate parameter rhohat from pairwise log-likelihoods
  rhohat_pw <- c()
  for(j in 1:p) {
    phi_pw <- phi[((j - 1)*d_rho + 1):(j*d_rho)]
    opt_pw <- nlminb(start = 0, objective = nloglik_beta, phi_st = phi_pw)
    rhohat_pw <- c(rhohat_pw, rho(opt_pw$par))
  }
  
  # i-th ubar by row - n x p
  # evaluated at rhohat
  Ubar <- t(apply(phi_mat_rho, 1, function(x) score_bar_rho(rho = rhohat_pw, 
                                                            phi_i = x)))
  lambda.val <- 56.3
  # lambda.val <- 16.95
  # lambda.val <- 2.94
  out.cell <- cd.alg(theta = rhohat_pw, lambda = lambda.val, scoremat = Ubar,
                     max.it = 30)
  w <- out.cell$what
  
  # selection
  nsel[i] <- sum(w != 0)
  
  rho_est <- rep(0, p)
  names(rho_est) <- names(rho_est)
  rho_est[w != 0] <- rhohat_pw[w != 0]
  
  sel <- which(rho_est != 0)
  
  rhotilde <- rep(0, p)
  rhotilde[sel] <- rhohat_pw[sel]
  
  rmse[i] <- sqrt(mean((rhohat - rhotilde)^2))
}

mean(rmse)
# [1] 0.2479353   ## 0.2938596    ### 0.1909254
mean(nsel)
# [1] 12        ## 6        ### 25

# pstarhat <- 12
# pstarhat <- 6
pstarhat <- 25

# Soft thresholding
rmse.thr.s <- rep(NA, ngroups)
nsel.thr.s <- rep(NA, ngroups)

for(i in 1:ngroups) {
  sub.i <- ss[!is.na(ss[, i]), i]
  y <- cell_mat[sub.i, ]
  n <- nrow(y)
  S <- cov(y)
  cv_soft <- threshold.cv(y, method = "soft", thresh.len = 5000, n.cv = 30,
                            norm = "F", seed = 321)
  thresh_soft <- cv_soft$threshold.grid
    
  sel_soft <- list(as.numeric(1000))
  nsel_soft <- rep(NA, length(thresh_soft))
  for(j in 1:length(thresh_soft)) {
    Sigma_soft <- soft.thresholding(S, threshold = thresh_soft[j])
    edges_soft <- as.vector(Sigma_soft[lower.tri(Sigma_soft)])
      
    sel_soft[[j]] <- which(edges_soft != 0)
    nsel_soft[j] <- length(which(edges_soft != 0))
  }
  ind.p <- which.min(abs(nsel_soft - pstarhat))
  Sigma_soft <- soft.thresholding(S, threshold = thresh_soft[ind.p])
  Sigma_soft <- cov2cor(Sigma_soft) 
  rhohat_soft <- as.vector(Sigma_soft[lower.tri(Sigma_soft)])
  
  nsel.thr.s[i] <- length(which(rhohat_soft != 0))
  
  names(rhohat_soft) <- names(rho_est)
  
  rmse.thr.s[i] <- sqrt(mean((rhohat - rhohat_soft)^2))
}

mean(rmse.thr.s)
# [1] 0.2921305   ## 0.3242812    ### 0.2270706
mean(nsel.thr.s)
# [1] 12          ## 6            ### 25.03333

# Bien & Tibshirani (2011) method
rmse.bt <- rep(NA, ngroups)
nsel.bt <- rep(NA, ngroups)

P <- matrix(1, d, d)            # equal penalty
diag(P) <- 0

# Set criterion for choice of lambda
pstar.tilde <- mean(nsel)

for(i in 1:ngroups) {
  sub.i <- ss[!is.na(ss[, i]), i]
  y <- cell_mat[sub.i, ]
  n <- nrow(y)
  S <- cov(y)
  
  ## lambdastar_bt <- 0.026
  lambdastar_bt <- 0.043
  # lambdastar_bt <- 0.0109
  Sigma_bt <- spcov(Sigma = S, S = S, lambda = lambdastar_bt * P,
                    step.size = 100)$Sigma
  Sigma_bt <- cov2cor(Sigma_bt)
  rhohat_bt <- as.vector(Sigma_bt[lower.tri(Sigma_bt)])
  nsel.bt[i] <- length(which(rhohat_bt != 0))
  
  names(rhohat_bt) <- names(rho_est)
  
  rmse.bt[i] <- sqrt(mean((rhohat - rhohat_bt)^2))
}

mean(rmse.bt)
# [1] 0.2922219  ## 0.34043     ### 0.2268005
mean(nsel.bt)
# [1] 12         ## 6           ### 25