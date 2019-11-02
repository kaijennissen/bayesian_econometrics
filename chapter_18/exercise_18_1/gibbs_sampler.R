library(Matrix)
library(svMisc)

rm(list = ls())
data <- read.csv("./chapter_18/exercise_18_1/USPCE_2015Q4.csv", header = FALSE)
y <- ts(data = data, start = c(1959, 1), freq = 4)
y <- c(diff(log(y))) * 400
y <- as.matrix(y)
colnames(y) <- NULL
TT <- dim(y)[1]

# priors
a0 <- 5
b0 <- 100

ny_sigma0 <- 3
S_sigma0 <- 2

ny_omega0 <- 3
S_omega0 <- 2 * 0.25^2


#
diags <- list(rep(1, TT), rep(-1, TT - 1))
H <- bandSparse(TT, k = c(0, -1), diag = diags, symm = FALSE)
HH <- t(H) %*% H # H is time invariant and therfore can be calculated outside the loop
HH <- as.matrix(HH)

nsim <- 1e4
nburn <- 1e3

store_tau <- matrix(0, nrow = nsim, ncol = TT)
store_tau0 <- matrix(0, nrow = nsim, ncol = 1)
store_sigma <- matrix(0, nrow = nsim, ncol = 1)
store_omega <- matrix(0, nrow = nsim, ncol = 1)

sigma2 <- 1
omega2 <- .1
tau0 <- 5
# tau <- y

bar_steps <- (nsim + nburn) / 100
for (ii in 1:(nsim + nburn)) {


  # draw from conditional for tau
  K_tau <- HH / omega2 + diag(TT) / sigma2
  B <- t(chol(K_tau))
  XX <- (HH * tau0) %*% matrix(1, nrow = TT) / c(omega2) + y / sigma2
  mu_tau <- solve(K_tau, XX)
  Z <- matrix(rnorm(n = TT, mean = 0, sd = 1))
  tau <- mu_tau + solve(t(B)) %*% Z

  # draw from conditional for sigma^2
  S_sigma_temp <- S_sigma0 + 0.5 * t(y - tau) %*% (y - tau)
  sigma2 <- 1 / rgamma(n = 1, shape = (ny_sigma0 + TT / 2), scale = 1 / S_sigma_temp)

  # draw from condition for omega^2
  S_omega_temp <- S_omega0 + 0.5 * t(tau - tau0) %*% HH %*% (tau - tau0)
  omega2 <- 1 / rgamma(n = 1, shape = ny_omega0 + TT / 2, scale = 1 / S_omega_temp)


  # draw from conditional for tau_0
  K_tau0 <- 1 / b0 + 1 / omega2
  invK_tau0 <- 1 / K_tau0
  tau0_hat <- invK_tau0 * (a0 / b0 + tau[1] / omega2)
  tau0 <- rnorm(n = 1, mean = tau0_hat, sd = sqrt(invK_tau0))

  progress(ii / bar_steps)
  if (ii == (nsim + nburn)) cat(": Done")

  if (ii > nburn) {
    nn <- ii - nburn
    store_tau[nn, ] <- tau[, , drop = T]
    store_tau0[nn, 1] <- tau0
    store_sigma[nn, 1] <- sigma2
    store_omega[nn, 1] <- omega2
  }
}

tau_post <- colMeans(store_tau)
sigma2_post <- colMeans(store_sigma)
omega2_post <- colMeans(store_omega)
tau0_post <- colMeans(store_tau0)


library(ggplot2)
library(lubridate)
library(tidyverse)
tt <- expand.grid(c(1959:2015), c(1:4))
xx <- lubridate::yq(paste(tt$Var1, tt$Var2, sep="-"))
tbl <- tibble(inflation = y, 'posterior mean' = colMeans(store_tau))
tbl <- mutate(tbl, t=sort(xx[-1]))


tbl %>% gather('series', 'value', -t) %>% 
ggplot(aes(x=t, y=value, lty = series, col=series)) +
  geom_line()
