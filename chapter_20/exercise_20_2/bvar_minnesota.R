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

us_macro <- read.csv("./chapter_20/exercise_20_2/US_macrodata.csv")
Yraw <- as.matrix(us_macro[!is.na(us_macro$INFLATION), 2:4])

tt <- nrow(Yraw) # number of observations
nn <- ncol(Yraw) # number of variables
intercept <- TRUE
pp <- 1 # lags, to be specified by the researcher
dd <- nn * (1 + nn * pp)

Ylag <- stats::lag(as.xts(ts(Yraw)), c(1:pp))

if (intercept) {
    X1 <- cbind(1, Ylag[(pp + 1):tt, ])
} else {
    X1 <- Ylag[(pp + 1):tt, ]
}

Y1 <- Yraw[(pp + 1):tt, ]

# -----------------------------
# OLS Estimates
# -----------------------------
Y <- Y1
X <- X1
TT <- dim(X)[1]
KK <- dim(X)[2]
A_OLS <- solve(t(X) %*% X) %*% (t(X) %*% Y)
a_OLS <- c(A_OLS)
SSE <- t(Y - X %*% A_OLS) %*% (Y - X %*% A_OLS)
Sigma_OLS <- SSE / (TT - KK)

y <- matrix(c(t(Yraw[-(1:pp), ])))

# -----------------------------
# Declare the prior values
# -----------------------------
mu0 <- matrix(0, nrow = dd, ncol = 1)
a1 <- 0.5
a2 <- 0.25
a3 <- 100
A1 <- matrix(a2, nrow = nn, ncol = nn)
diag(A1) <- a1
A2 <- matrix(rep(1, pp), ncol = pp) %x% A1
A <- c(cbind(a3, A2))

# -----------------------------
# Calculate Variance Matrix
# -----------------------------
DSigma_OLS <- diag(diag(Sigma_OLS))
Sig1 <- bdiag(DSigma_OLS, solve(DSigma_OLS) %x% diag((nn * pp))) # %x% DSigma_OLS
invSig1 <- bdiag(diag(nn), solve(DSigma_OLS) %x% diag((nn * pp)))
R <- bdiag(diag(nn), diag((1:pp)^-2) %x% diag(nn^2))
V0 <- diag(A) %*% Sig1 %*% invSig1 %*% R
invV0 <- solve(V0)

#---------------------------------
# Calculate the Posterior
#---------------------------------
for (j in (pp + 1):tt) {
    X_temp <- diag(nn) %x% cbind(1, t(c(t(Yraw[(j - 1):(j - pp), , drop = FALSE]))))
    if (j == (pp + 1)) {
        X <- X_temp
    } else {
        X <- rbind(X, X_temp)
    }
}

XiSig <- t(X) %*% (diag((tt - pp)) %x% solve(Sigma_OLS))
K_beta <- XiSig %*% X + invV0
beta_hat <- solve(K_beta) %*% (XiSig %*% y)
b <- matrix(beta_hat[1:nn])
A1 <- matrix(beta_hat[(nn + 1):(nn * (1 + nn))], ncol = nn)

# ---------------------------------------------------------------------------------------
# # Y(t) =  b + A1 * Y(t-1) +  A2 * Y(t-2) + ... + Ap * Y(t-p)
# ---------------------------------------------------------------------------------------
