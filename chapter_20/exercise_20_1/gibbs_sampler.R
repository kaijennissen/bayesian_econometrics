# -----------------------------------------------------------------------------
# R Code for Exercise 20.1
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
# -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)

us_macro <- read.csv("./chapter_20/exercise_20_1/US_macrodata.csv")
us_macro <- as.matrix(us_macro[!is.na(us_macro$INFLATION), 2:4])

tt <- nrow(us_macro) # number of observations
nn <- ncol(us_macro) # number of variables
pp <- 2 # lags; to be specified by the researcher
dd <- nn * (1 + nn * pp)

y <- matrix(c(t(us_macro[-(1:(pp)), ])))
# assertthat::are_equal(dim(y)[1], (tt-pp)*(nn))
# assertthat::are_equal(dim(y)[2],1)

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


niter <- 2e4
burn <- 1e3

#---------------------------------
# Set lengths of parameter vectors
#---------------------------------
store_beta <- matrix(NA, nrow = dd, ncol = niter)
store_Sigma <- array(NA, dim = c(nn, nn, niter))

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
  }
}

round(apply(store_beta, MARGIN = c(1), FUN = mean), 4)
round(apply(store_Sigma, MARGIN = c(1, 2), FUN = mean), 4)
