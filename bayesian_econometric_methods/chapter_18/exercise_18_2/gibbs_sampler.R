#------------------------------------------------------------------------------
# R Code for Exercise 18.2
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
#------------------------------------------------------------------------------
# run time: 96.6 sec

library(Matrix)
library(dplyr)
library(ggplot2)
library(RcppZiggurat)
library(Rcpp)

rm(list = ls())
now <- Sys.time()

set.seed(123)

nsim <- 10000
nburn <- 1000
total_runs <- nsim + nburn

data <-
    read.csv("./bayesian_econometric_methods/chapter_18/exercise_18_2/usgdp.csv", header = FALSE)
y <- ts(data = data,
        start = c(1959, 1),
        freq = 4)
y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL

TT <- dim(y)[1]
p <- 2



# priors #---------------------------------------------------------------------
# phi
phi_0 <- matrix(c(1.3, -0.7))
V_phi <- diag(2)
invV_phi <- solve(V_phi)

# gamma = (tau(0), tau(-1))'
gamma_0 <- matrix(c(750, 750))
V_gamma <- 100 * diag(2)
invV_gamma <- solve(V_gamma)

# sigma2_tau
b_sigma2_tau <- 0.01

# sigma2_c
ny_sigma2_c <- 3
S_sigma2_c <- 2


# initial values #-------------------------------------------------------------

f_tau <- function(x, TT, del_tau) {
    -TT / 2 * log(x) - sum(diff(del_tau)[-TT] ** 2) / (2 * x)
}


phi <- Matrix(c(1.34, -0.7))
sigma2_c <- .5
sigma2_tau <- .001
gamma <- matrix(c(y[1], y[1])) # !! gamma = (tau(0), tau(-1))'

# H_2
diags_h2 <- list(rep(1, TT), rep(-2, TT - 1), rep(1, TT - 2))
H_2 <-
    bandSparse(TT,
               k = c(0, -1, -2),
               diag = diags_h2,
               symm = FALSE)


invH_2 <- solve(H_2, diag(TT))
HH_2 <- crossprod(H_2)


# H_phi
diags_phi <-
    list(rep(1, TT), rep(-phi[1], TT - 1), rep(-phi[2], TT - 2))
H_phi <-
    bandSparse(TT,
               k = c(0, -1, -2),
               diag = diags_phi,
               symm = FALSE)
#HH_phi <- t(H_phi) %*% H_phi
HH_phi <- crossprod(H_phi)

# X_gamma
X_gamma <- cbind(c(2:(TT + 1)), -c(1:TT))

store_tau <- matrix(0, nrow = TT, ncol = nsim)
store_theta <- matrix(0, nrow = 6, ncol = nsim)

n_grid <- 500
count_phi <- 0

alpha_tau_tilde <- matrix(0, nrow = TT)


for (ii in 1:total_runs) {
    # draw from conditional for tau #------------------------------------------
    # gamma = (tau(0), tau(-1))'
    alpha_tau_tilde[1, 1] <- 2 * gamma[1] - gamma[2]
    alpha_tau_tilde[2, 1] <- -gamma[1]
    alpha_tau <- invH_2 %*% alpha_tau_tilde
    K_tau <- HH_phi / sigma2_c + HH_2 / sigma2_tau
    XX_tau <-
        (HH_phi %*% y) / sigma2_c + (HH_2 %*% alpha_tau) / sigma2_tau
    L <- chol(K_tau)
    tau_hat <- Matrix(backsolve(L, forwardsolve(L, XX_tau, upper.tri = FALSE, transpose = FALSE), upper.tri = TRUE, transpose = FALSE))
    tau <- tau_hat + solve(L, zrnorm(n = TT))
    tau <- tau@x


    
    # draw from conditional for phi #------------------------------------------
    cc <- y - tau
    
    X_phi <- cbind(c(0, cc[1:(TT - 1)]), c(0, 0, cc[1:(TT - 2)]))
    XX_phi <- crossprod(X_phi)
    K_phi <- invV_phi + XX_phi / sigma2_c
    L_phi <- chol(K_phi) #
    XX <- invV_phi %*% phi_0 + crossprod(X_phi, cc) / sigma2_c
    phi_hat <- solve(K_phi, XX)
    #Z_phi <- matrix(rnorm(n = p, mean = 0, sd = 1))
    # Z_phi <- matrix(rziggurat(n = p))
    Z_phi <- zrnorm(n = p)
    phic <- phi_hat + solve(L_phi) %*% Z_phi
    # stationarity region for AR(2)
    # 1) phi_1 + phi_2 < 1
    # 2) phi_2 - phi_1 < 1
    # 3) phi_2 > -1
    if (sum(phic) < .99 &&
        phic[2] - phic[1] < .99 && phic[2] > -.99) {
        phi <- phic
        # calculate H_phi in every iteration
        diags_phi <-
            list(rep(1, TT), rep(-phi[1], TT - 1), rep(-phi[2], TT - 2))
        H_phi <-
            bandSparse(TT,
                       k = c(0, -1, -2),
                       diag = diags_phi,
                       symm = FALSE)
        HH_phi <- crossprod(H_phi)
        count_phi <- count_phi + 1
    }
    
    
    # draw from condition for sigma^2_c #--------------------------------------
    S_sigma2_c_temp <- S_sigma2_c + 0.5 * t(cc) %*% (HH_phi %*% cc)
    sigma2_c <-
        1 / rgamma(
            n = 1,
            shape = ny_sigma2_c + TT / 2,
            scale = 1 / S_sigma2_c_temp@x
        )
    
    
    # draw from conditional for sigma^2_tau #---------------------------------------------------------
    # gamma = (tau(0), tau(-1))'
    del_tau <-
        c(gamma[1], tau[1:TT]) - c(gamma[2], gamma[1], tau[1:(TT - 1)])
    sigma2_tau_grid <-
        seq(
            from = runif(1) / 1000,
            to = b_sigma2_tau - runif(1) / 1000,
            length.out = n_grid
        )
    lp_sigtau2 <-
        f_tau(x = sigma2_tau_grid,
              TT = TT,
              del_tau = del_tau)
    p_sigtau2 <- exp(lp_sigtau2 - max(lp_sigtau2))
    p_sigtau2 <- p_sigtau2 / sum(p_sigtau2)
    cdf_sigtau2 <- cumsum(p_sigtau2)
    sigma2_tau <-
        sigma2_tau_grid[which(runif(1) < cdf_sigtau2, TRUE)[1]]
    
    
    # draw from conditional for gamma #----------------------------------------
    K_gamma <-
        invV_gamma + t(X_gamma) %*% (HH_2 %*% X_gamma) / sigma2_tau
    L_gamma <- chol(K_gamma)
    XX <-
        invV_gamma %*% gamma_0 + t(X_gamma) %*% (HH_2 %*% tau) / sigma2_tau
    gamma_hat <- solve(K_gamma, XX)
    #Z_gamma <- rnorm(n = p, mean = 0, sd = 1)
    Z_gamma <- zrnorm(n = p)
    gamma <- gamma_hat + solve(L_gamma) %*% Z_gamma
    
    
    # store #------------------------------------------------------------------
    if (ii > nburn) {
        nn <- ii - nburn
        store_tau[, nn] <- tau
        store_theta[, nn] <-
            c(c(phi), sigma2_c, sigma2_tau, gamma@x)
    }
}
print(Sys.time() - now)


# Results #--------------------------------------------------------------------
tau_hat <- rowMeans(store_tau)
theta_hat <- rowMeans(store_theta)
print(theta_hat)
# 1.313774069  -0.380283872   0.811751869   0.002302012 747.395159671 746.471692513

cc <- ts(y - tau_hat, start = c(1949, 1), frequency = 4)

timetk::tk_tbl(data = cc, rename_index = "time") %>%
    rename(c = "Series 1") %>%
    ggplot2::ggplot(aes(x = time, y = c)) +
    geom_line()

