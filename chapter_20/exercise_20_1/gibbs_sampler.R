# -----------------------------------------------------------------------------
# R Code for Exercise 20.1
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
# -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)

# -----------------------------
# Impulse Resonse
# -----------------------------
construct_IR <- function(beta, Sig, n_hz, shock) {
  nn <- dim(Sig)[1]
  pp <- (dim(beta)[1] / nn - 1) / nn
  CSig <- t(chol(Sig))
  tmpZ1 <- matrix(0, nrow = pp, ncol = nn)
  tmpZ <- matrix(0, nrow = pp, ncol = nn)
  Yt1 <- CSig %*% shock
  Yt <- matrix(0, nrow = nn)
  yIR <- matrix(0, nrow = n_hz, ncol = nn)
  yIR[1, ] <- t(Yt1)
  for (t in 2:n_hz) {
    # update the regressors
    len <- dim(tmpZ)[1]
    tmpZ <- rbind(t(Yt), tmpZ[1:(len - 1), ])
    tmpZ1 <- rbind(t(Yt1), tmpZ1[1:(len - 1), ])
    # evolution of variables if a shock hits
    e <- CSig %*% matrix(rnorm(nn, 1))
    Z1 <- matrix(c(t(tmpZ1)), nrow = 1, ncol = nn * pp)
    Xt1 <- diag(nn) %x% cbind(1, Z1)
    Yt1 <- Xt1 %*% beta + e
    # evolution of variables if no shocks hit

    Z <- matrix(c(t(tmpZ)), nrow = 1, ncol = nn * pp)
    Xt <- diag(nn) %x% cbind(1, Z)
    Yt <- Xt %*% beta + e
    # the IR is the difference of the two scenarios
    yIR[t, ] <- c(t(Yt1 - Yt))
  }
  return(yIR)
}

us_macro <- read.csv("./chapter_20/exercise_20_1/US_macrodata.csv")
us_macro <- as.matrix(us_macro[!is.na(us_macro$INFLATION), 2:4])

tt <- nrow(us_macro) # number of observations
nn <- ncol(us_macro) # number of variables
pp <- 2 # lags; to be specified by the researcher
dd <- nn * (1 + nn * pp)
n_hz <- 20

y <- matrix(c(t(us_macro[-(1:(pp)), ])))

# -----------------------------
# Declare the prior values
# -----------------------------
# priors for beta
mu0 <- matrix(0, nrow = dd, ncol = 1)
Sigma0 <- diag(dd)
diag(Sigma0) <- rep(10, dd)
invSigma0 <- solve(Sigma0)

# priors for sigma
v0 <- nn + 3
S0 <- diag(nn)
diag(S0) <- rep(1, nn)


niter <- 2e3
burn <- 1e3

#---------------------------------
# Set lengths of parameter vectors
#---------------------------------
store_beta <- matrix(NA, nrow = dd, ncol = niter)
store_Sigma <- array(NA, dim = c(nn, nn, niter))
store_yIR <- matrix(0, nrow = n_hz, ncol = nn)

#----------------------
# Set initial values
#----------------------
for (j in (pp + 1):tt) {
  X_temp <- diag(nn) %x% cbind(1, t(c(t(us_macro[(j - 1):(j - pp), , drop = FALSE]))))
  if (j == (pp + 1)) {
    X <- X_temp
  } else {
    X <- rbind(X, X_temp)
  }
}

# initialize the chain
beta <- solve(t(X) %*% X) %*% (t(X) %*% y)
e <- matrix(y - X %*% beta, nrow = nn, ncol = (tt - pp))
Sigma <- e %*% t(e) / tt


#----------------------
# Begin the Gibbs Sampler
#----------------------
for (i in 1:(niter + burn)) {
  # Do conditional for beta
  XiSig <- t(X) %*% (diag((tt - pp)) %x% solve(Sigma))
  Dtemp <- solve(XiSig %*% X + invSigma0)
  dtemp <- (XiSig %*% y + invSigma0 %*% mu0)
  LL <- chol(Dtemp) # chol returns the upper triangular
  beta <- Dtemp %*% dtemp + t(LL) %*% matrix(rnorm(dd))

  # Do conditional for Sigma
  Dtemp <- matrix(0, nrow = nn, ncol = nn)
  for (j in (pp + 1):tt) {
    X_temp <- diag(nn) %x% cbind(1, t(c(t(us_macro[(j - 1):(j - pp), , drop = FALSE]))))
    y_temp <- t(us_macro[j, , drop = FALSE])
    eps_temp <- y_temp - X_temp %*% beta
    Dtemp <- Dtemp + eps_temp %*% t(eps_temp)
  }
  Dtemp <- Dtemp + S0
  dtemp <- v0 + (tt - pp)

  Sigma <- solve(rWishart(1, Sigma = solve(Dtemp), df = dtemp)[, , 1])
  if (i %% 1000 == 0) {
    print(paste("loops... ", i))
  }

  if (i > burn) {
    nsim <- i - burn
    store_beta[, nsim] <- beta
    store_Sigma[, , nsim] <- Sigma

    # calculate impulse response
    CSig <- t(chol(Sigma))
    # 100 basis pts rather than 1 std. dev.
    shock <- matrix(c(0, 0, 1)) / CSig[nn, nn]
    yIR <- construct_IR(beta = beta, Sig = Sigma, n_hz = n_hz, shock = shock)

    store_yIR <- store_yIR + yIR
  }
}

beta_hat <- round(apply(store_beta, MARGIN = c(1), FUN = mean), 4)
sigma_hat <- round(apply(store_Sigma, MARGIN = c(1, 2), FUN = mean), 4)
yIR_hat <- store_yIR / niter

par(mfrow = c(3, 1))
plot.ts(yIR_hat[, 1], ylab = "Inflation", main = "Impulse Response")
plot.ts(yIR_hat[, 2], ylab = "Unemployment Rate")
plot.ts(yIR_hat[, 3], ylab = "Fed Funds Rate")
