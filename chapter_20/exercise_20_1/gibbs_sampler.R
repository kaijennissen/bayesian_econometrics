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
mu_beta <- matrix(0, nrow = dd, ncol = 1)
# assertthat::are_equal(dim(mu_beta)[1],dd)
# assertthat::are_equal(dim(mu_beta)[2],1)
sigma_beta <- diag(dd)
diag(sigma_beta) <- rep(10, dd)
# assertthat::are_equal(dim(mu_beta)[1],dd)
# assertthat::are_equal(dim(mu_beta)[2],dd)

v0 <- nn+3
S0 <- diag(nn)
diag(S0) <- rep(1, nn)
# assertthat::are_equal(dim(S0)[1],nn)
# assertthat::are_equal(dim(S0)[2],nn)

niter <- 2e4
burn <- 1e3

#---------------------------------
# Set lengths of parameter vectors
#---------------------------------
beta <- matrix(NA, nrow = dd, ncol = niter)
Sigma <- array(NA, dim = c(nn, nn, niter))

#----------------------
# Set initial conditions
#----------------------
beta[, 1] <- matrix(rnorm(dd), nrow = dd, ncol = 1)
Sigma[, , 1] <- matrix(rgamma(n = nn * nn, shape = 2, scale = 2), nrow = nn, ncol = nn) + 111 * diag(nn)
# chol(Sigma[,,1])

#----------------------
# Begin the Gibbs Sampler
#----------------------
# i=2
for (i in 2:niter) {
  # Do conditional for beta
  for (j in (pp + 1):tt) {
    X_temp <- diag(nn) %x% cbind(1, t(c(t(us_macro[(j - 1):(j - pp), , drop = FALSE]))))
    if (j == (pp + 1)) {
      X <- X_temp
    } else {
      X <- rbind(X, X_temp)
    }
  }
  # assertthat::are_equal(dim(X)[1], (tt-pp)*nn)
  # assertthat::are_equal(dim(X)[2], dd)

  dtemp <- t(X) %*% (diag((tt - pp)) %x% solve(Sigma[, , i - 1])) %*% y + solve(sigma_beta) %*% mu_beta
  Dtemp <- solve(t(X) %*% (diag((tt - pp)) %x% solve(Sigma[, , i - 1])) %*% X + solve(sigma_beta))
  H <- chol(Dtemp) # chol returns the upper triangular, therefore use tranposed in next step
  beta[, i] <- Dtemp %*% dtemp + t(H) %*% matrix(rnorm(dd))

  # Do conditional for Sigma
  Dtemp <- matrix(0, nrow = nn, ncol = nn)
  for (j in (pp + 1):tt) {
    X_temp <- diag(nn) %x% cbind(1, t(c(t(us_macro[(j - 1):(j - pp), , drop = FALSE]))))
    y_temp <- t(us_macro[j, , drop = FALSE])
    eps_temp <- y_temp - X_temp %*% beta[, i]
    Dtemp <- Dtemp + eps_temp %*% t(eps_temp)
  }
  Dtemp <- Dtemp + S0
  dtemp <- v0 + (tt - pp)

  #  Sigma[, , i] <- rWishart(n = 1, df = dtemp*Dtemp, Sigma = Dtemp)
  Sigma[, , i] <- solve(rWishart(1, Sigma = Dtemp, df = dtemp)[, , 1]) #  i=i+1
  if (i %% 1000 == 0) {
    print(paste("loops... ", i))
  }
  
  # # compute impulse-responses
  # CSig <-  t(chol(Sig[,,1]))
  # #100 basis pts rather than 1 std. dev.
  # shock <- matrix(c(0, 0, 1))%*%solve(CSig[nn,nn])
  # yIR <- construct_IR(beta,Sig,n_hz,shock)
  # store_yIR <-  store_yIR + yIR
}

beta_post <- beta[, -(1:burn)]
Sigma_post <- Sigma[, , -(1:burn)]

apply(beta_post, MARGIN = c(1), FUN = mean)
apply(Sigma_post, MARGIN = c(1, 2), FUN = mean)

