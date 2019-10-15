# -----------------------------------------------------------------------------
# R Code for Exercise 20.2
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
# -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)
library(xts)
library(Matrix)
# -----------------------------
# Impulse Resonse
# -----------------------------

us_macro <- read.csv("./chapter_20/exercise_20_2/US_macrodata.csv")
Yraw <- as.matrix(us_macro[!is.na(us_macro$INFLATION), 2:4])

tt <- nrow(Yraw) # number of observations
nn <- ncol(Yraw) # number of variables
pp <- 1 # lags; to be specified by the researcher
intercept <- FALSE
dd <- nn * (1 + nn * pp)
n_hz <- 20

# 
Ylag <- stats::lag(as.xts(ts(Yraw)), c(1:pp))

if (intercept) {
  X1 <- cbind(1, Ylag[(pp + 1):tt, ])
} else {
  X1 <- Ylag[(pp + 1):tt, ]
}

Y1 <-  Yraw[(pp+1):tt,]

# -----------------------------
# OLS Estimates
# -----------------------------
Y <- Y1
X <- X1
TT <- dim(X)[1]
KK <- dim(X)[2]
A_OLS <-  solve(t(X)%*%X)%*%(t(X)%*%Y)
#a_OLS <-  A_OLS(:);         #% This is the vector of coefficients, i.e. it holds vec(A_OLS)
SSE <-  t(Y - X%*%A_OLS)%*%(Y - X%*%A_OLS)
Sigma_OLS <-  SSE/(TT-KK)



y <- matrix(c(t(Yraw[-(1:pp), ])))

# -----------------------------
# Declare the prior values
# -----------------------------
# priors for beta
mu0 <- matrix(0, nrow = dd, ncol = 1)
a1 <- 0.5
a2 <- 0.25
a3 <- 100
A1 <- matrix(a2, nrow = nn, ncol = nn)
diag(A1) <- a1
A2 <- matrix(rep(1, pp), ncol = pp) %x% A1
A <- c(cbind(a3, A2))

# Sigma0 <- diag(nn)
# diag(Sigma0) <- rep(10, dd)
# invSigma0 <- solve(Sigma0)

DSigma_OLS <- diag(diag(Sigma_OLS))
Sig1 <- bdiag(DSigma_OLS, solve(DSigma_OLS) %x% diag((nn*pp))) #%x% DSigma_OLS
invSig1 <- bdiag(diag(nn), solve(DSigma_OLS) %x% diag((nn*pp)))
R <- bdiag(diag(nn), diag((1:pp)^-2)%x%diag(nn^2))
V0 <- diag(A) %*% Sig1 %*% invSig1 %*% R
invV0 <- solve(V0)

# priors for sigma
# v0 <- nn + 3
# S0 <- diag(nn)
# diag(S0) <- rep(1, nn)


niter <- 2e4
burn <- 1e3


mlag2 <- function(X, p) {
  X <- as.matrix(X)
  # we need to bind horizontally p blocks
  Xlag <- matrix(nrow = nrow(X), ncol = 0)  # create empty matrix with correct number of raws
  for (i in 1:p) {
    Xlag <- cbind(Xlag, makeblock(X, i, p))  # bind blocks horizontally
  }
  return(Xlag)
}


makeblock <- function(X, i, p) {
  Xblock <- X[(p + 1 - i):(nrow(X) - i), ]  # get useful lines of X
  Xblock <- as.matrix(Xblock)  # assure X is a matrix, not a vector
  Xblock <- rbind(matrix(0, nrow = p, ncol = ncol(X)), Xblock)  # append p zero lines at top
  return(Xblock)
}


#---------------------------------
# Set lengths of parameter vectors
#---------------------------------
store_beta <- matrix(NA, nrow = dd, ncol = niter)
#store_Sigma <- array(NA, dim = c(nn, nn, niter))
store_yIR <- matrix(0, nrow = n_hz, ncol = nn)

#----------------------
# Set initial values
#----------------------
for (j in (pp + 1):tt) {
  X_temp <- diag(nn) %x% cbind(1, t(c(t(Yraw[(j - 1):(j - pp), , drop = FALSE]))))
  if (j == (pp + 1)) {
    X <- X_temp
  } else {
    X <- rbind(X, X_temp)
  }
}

#---------------------------------
# Posterior
#---------------------------------
XiSig <-  t(X)%*%(diag((tt-pp))%x%solve(Sigma_OLS))
K_beta <-XiSig%*%X+invV0
beta_hat <- solve(K_beta)%*%(XiSig%*%y)
#b <- matrix(beta_hat[1:nn])
A1 <- matrix(beta_hat[4:12], ncol=nn)

# initialize the chain
beta <- matrix(0, nrow=nn*(1+nn*pp))#solve(t(X) %*% X) %*% (t(X) %*% y)
#e <- matrix(y - X %*% beta, nrow = nn, ncol = (tt - pp))
Sigma <- Sigma_OLS
#invSigma <- solve(Sigma_OLS)



#----------------------
# Begin the Gibbs Sampler
#----------------------
for (i in 1:(niter + burn)) {
  # Do conditional for beta
  XiSig <- t(X) %*% (diag((tt - pp)) %x% solve(Sigma))
  Dtemp <- solve(XiSig %*% X + invV0)
  dtemp <- (XiSig %*% y + invV0 %*% mu0)
  LL <- chol(as.matrix(Dtemp)) # chol returns the upper triangular
  beta <- Dtemp %*% dtemp + t(LL) %*% matrix(rnorm(dd))

  # Do conditional for Sigma
  # Dtemp <- matrix(0, nrow = nn, ncol = nn)
  # for (j in (pp + 1):tt) {
  #   X_temp <- diag(nn) %x% cbind(1, t(c(t(us_macro[(j - 1):(j - pp), , drop = FALSE]))))
  #   y_temp <- t(us_macro[j, , drop = FALSE])
  #   eps_temp <- y_temp - X_temp %*% beta
  #   Dtemp <- Dtemp + eps_temp %*% t(eps_temp)
  # }
  # Dtemp <- Dtemp + S0
  # #dtemp <- v0 + (tt - pp)
  #Sigma <- solve(rWishart(1, Sigma = solve(Dtemp), df = dtemp)[, , 1])
  
  if (i %% 1000 == 0) {
    print(paste("loops... ", i))
  }

  if (i > burn) {
    nsim <- i - burn
    store_beta[, nsim] <- as.matrix(beta)
    #store_Sigma[, , nsim] <- Sigma

    # calculate impulse response
    CSig <- t(chol(Sigma))
    # 100 basis pts rather than 1 std. dev.
    shock <- matrix(c(0, 0, 1)) / CSig[nn, nn]
    yIR <- construct_IR(beta = beta, Sig = Sigma, n_hz = n_hz, shock = shock)

    store_yIR <- store_yIR + yIR
  }
}

beta_hat <- round(apply(store_beta, MARGIN = c(1), FUN = mean), 4)
b <- matrix(beta_hat[1:nn])
A1 <- matrix(beta_hat[4:12], ncol=nn)

# -----------------------------
# # Y(t) =  b + A(1) * Y(t-1)
# -----------------------------
