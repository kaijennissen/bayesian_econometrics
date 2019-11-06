# -----------------------------------------------------------------------------
# R Code for Exercise 18.2
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
# -----------------------------------------------------------------------------
library(Matrix)
library(progress)
library(tictoc)
# library(Rfast)

rm(list = ls())


#data <- read.csv("./chapter_18/exercise_18_2/USPCE_2015Q4.csv", header = FALSE)
data <- read.csv("./chapter_18/exercise_18_2/usgdp.csv", header = FALSE)
y <- ts(data = data, start = c(1959, 1), freq = 4)
y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL
TT <- dim(y)[1]
p <- 2

#-------------------
# priors
#-------------------
# phi
phi_0 <- matrix(c(1.3, -0.7))
V_phi <- diag(2)
invV_phi <- Rfast::spdinv(V_phi)

# gamma
gamma_0 <- matrix(c(750, 750))
V_gamma <- 100 * diag(2)
invV_gamma <- Rfast::spdinv(V_gamma)

# sigma2_tau
b_sigma2_tau <- 0.01

# sigma2_c
ny_sigma2_c <- 3
S_sigma2_c <- 2

#-------------------
# initla values
#-------------------
phi <- matrix(c(1.34, -0.7))
sigma2_c <- .5
sigma2_tau <- .001
tau0 <- matrix(c(y[1], y[1]))

# H_2
diags <- list(rep(1, TT), rep(-2, TT - 1), rep(1, TT - 2))
H_2 <- as.matrix(Matrix::bandSparse(TT, k = c(0, -1, -2), diag = diags, symm = FALSE))
invH_2 <- spdinv(H_2)
HH_2 <- Rfast::mat.mult(Rfast::transpose(H_2), H_2) #t(H_2) %*% H_2


# H_phi
diags_phi <- list(rep(1, TT), rep(-phi[1], TT - 1), rep(-phi[2], TT - 2))
H_phi <- as.matrix(Matrix::bandSparse(TT, k = c(0, -1, -2), diag = diags_phi, symm = FALSE))
HH_phi<- Rfast::mat.mult(t(H_phi), H_phi)#  crossprod(H_phi, H_phi) #t(H_phi) %*% H_phi

# X_gamma
X_gamma <- cbind(c(2:(TT + 1)), -c(1:TT))


nsim <- 1e1
nburn <- 1e1

store_tau <- matrix(0, nrow = nsim, ncol = TT)
store_gamma <- matrix(0, nrow = nsim, ncol = 2)
store_sigma2_tau <- matrix(0, nrow = nsim, ncol = 1)
store_sigma2_c <- matrix(0, nrow = nsim, ncol = 1)
store_phi <- matrix(0, nrow = nsim, ncol = 2)

n_grid <- 500
count_phi <- 0

step_size <- (nsim + nburn) / 100
pb <- progress_bar$new(
  format = "[:bar] :percent in :elapsed",
  total = nsim+nburn, clear = FALSE, width = 60
)


for (ii in 1:(nsim + nburn)) {
    if(ii==1) cat("\014")
    pb$tick()
    Sys.sleep(1 / 100)

  
  # draw from conditional for tau
  alpha_tau_tilde <- matrix(c(2 * tau0[2] - tau0[1], -tau0[2], rep(0, TT - 2)))
  alpha_tau <- Rfast::mat.mult(invH_2, alpha_tau_tilde)#invH_2 %*% alpha_tau_tilde
  K_tau <- HH_phi / sigma2_c + HH_2 / sigma2_tau
  L_tau <- Rfast::cholesky(K_tau) # returns upper
  XX_tau <-Rfast:: mat.mult(HH_phi, y) / sigma2_c + Rfast::mat.mult(HH_2, alpha_tau) / sigma2_tau
  tau_hat <- solve(K_tau, XX_tau)
  Z_tau <- matrix(Rfast::Rnorm(n = TT, m = 0, s = 1))
  tau <- tau_hat + Rfast::mat.mult(chol2inv(K_tau), Z_tau)


  # draw from conditional for phi
  c <- y - tau
  X_phi <- cbind(c(0, c[1:(TT - 1)]), c(0, 0, c[1:(TT - 2)]))
  XX_phi <- Rfast::mat.mult(Rfast::transpose(X_phi), X_phi)
  K_phi <- invV_phi + XX_phi / sigma2_c
  L_phi <- Rfast::cholesky(K_phi) #
  XX <- Rfast::mat.mult(invV_phi, phi_0) + Rfast::mat.mult(Rfast::transpose(X_phi), c) / sigma2_c
  phi_hat <- solve(K_phi, XX)
  Z_phi <- matrix(Rfast::Rnorm(n = p, m = 0, s = 1))
  phic <- phi_hat + Rfast::mat.mult(chol2inv(K_phi), Z_phi)
  # stationarity region for AR(2)
  # 1) phi_1 + phi_2 < 1
  # 2)
  if (sum(phic) < .99 && phic[2] - phic[1] < .99 && phic[2] > -.99) {
    phi <- phic
    # calculate H_phi in every iteration
    diags_phi <- list(rep(1, TT), rep(-phi[1], TT - 1), rep(-phi[2], TT - 2))
    H_phi <-  as.matrix(bandSparse(TT, k = c(0, -1, -2), diag = diags_phi, symm = FALSE))
    HH_phi <-Rfast::mat.mult(Rfast::transpose(H_phi),  H_phi)
    count_phi <- count_phi + 1
  }


  # draw from condition for sigma^2_c
  S_sigma2_c_temp <- S_sigma2_c + 0.5 * Rfast::mat.mult(Rfast::mat.mult(Rfast::transpose(y - tau), HH_phi), (y - tau))
  sigma2_c <- 1 / rgamma(n = 1, shape = ny_sigma2_c + TT / 2, scale = 1 / S_sigma2_c_temp)


  # draw from conditional for sigma^2_tau
  del_tau <- c(tau0[1], tau[1:TT]) - c(tau0[2], tau0[1], tau[1:(TT - 1)])
  f_tau <- function(x) {
    -TT / 2 * log(x) - sum(diff(del_tau)[-TT]**2) / (2 * x)
  }
  # sum((del_tau[2:TT]-del_tau[1:(TT-1)])**2)
  sigma2_tau_grid <- seq(from = runif(1) / 1000, to = b_sigma2_tau - runif(1) / 1000, length.out = n_grid)
  lp_sigtau2 <- f_tau(sigma2_tau_grid)
  p_sigtau2 <- exp(lp_sigtau2 - max(lp_sigtau2))
  p_sigtau2 <- p_sigtau2 / sum(p_sigtau2)
  cdf_sigtau2 <- cumsum(p_sigtau2)
  sigma2_tau <- sigma2_tau_grid[which(runif(1) < cdf_sigtau2, TRUE)[1]]


  # draw from conditional for gamma
  K_gamma <- invV_gamma + Rfast::mat.mult(Rfast::transpose(X_gamma),  Rfast::mat.mult(HH_2, X_gamma) )/ sigma2_tau
  L_gamma <- chol(K_gamma)
  XX <- Rfast::mat.mult(invV_gamma, gamma_0) +  Rfast::mat.mult(Rfast::transpose(X_gamma),  Rfast::mat.mult(HH_2, tau)) / sigma2_tau
  gamma_hat <- solve(K_gamma, XX)
  Z_gamma <- matrix(Rfast::Rnorm(n = p, m = 0, s = 1))
  gamma <- gamma_hat + Rfast::mat.mult(solve(L_gamma), Z_gamma)

  if (ii > nburn) {
    nn <- ii - nburn
    store_tau[nn, ] <- tau[,,drop=TRUE]
    store_phi[nn, ] <- phi[,,drop=TRUE]
    store_gamma[nn, ] <- gamma[,,drop=TRUE]
    store_sigma2_tau[nn] <- sigma2_tau
    store_sigma2_c[nn] <- sigma2_c
  }
}

tau_post <- colMeans(store_tau)
phi_post <- colMeans(store_phi)
gamma_post <- colMeans(store_gamma)
sigma2_tau_post <- colMeans(store_sigma2_tau)
sigma2_c_post <- colMeans(store_sigma2_c)




## Rjags

library(rjags)
library(coda)
























