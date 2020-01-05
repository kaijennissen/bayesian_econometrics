# -----------------------------------------------------------------------------
# R Code for Exercise 16.7
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
# -----------------------------------------------------------------------------

rm(list = ls())
set.seed(123)

# -----------------------------
# Simulate data
# -----------------------------
Sigma_sim <- diag(5)
Sigma_sim[lower.tri(Sigma_sim)] <- c(.4, .4, 0, .6, .7, 0, .3, 0, .3, 0)
Sigma_sim[upper.tri(Sigma_sim)] <- t(Sigma_sim)[upper.tri(Sigma_sim)]

L <- t(chol(Sigma_sim))
X_sim <- cbind(1, t(L %*% matrix(rnorm(5000), nrow = 5)))
beta_sim <- c(2, .25, -.4, .6, 0, 0)

y_sim <- X_sim %*% beta_sim + 0.2 * rnorm(1000)

X <- X_sim
y <- y_sim


# -----------------------------
# SSVS
# -----------------------------
nn <- nrow(y_sim)
kk <- 5
niter <- 1e5 # declare the number of iterations
burn <- 2e2 # declare the length of the burn-in

# -----------------------------
# Declare the prior values
# -----------------------------
# prior for gamma_i
p <- 0.5
gammas <- matrix(0, nrow = 5)
# note that for a small number of iterations the posterior probability
# depend heavily on the initial values

# prior for beta
mu_sigeps <- .15
tau2 <- 1e-8
c2 <- 9 / tau2
V0 <- 10^2

# prior for sigma^2
a <- 3
b <- 1 / (2 * mu_sigeps)
sigsq <- 1


#---------------------------------
# Store samples
#---------------------------------
store_beta <- matrix(0, nrow = kk + 1, ncol = niter)
store_sigma <- matrix(0, nrow = 1, ncol = niter)
store_gamma <- matrix(0, nrow = kk, ncol = niter)

#----------------------
# Begin the Gibbs Sampler
#----------------------
for (i in 1:(niter + burn)) {
    # i=i+1
    # Conditional for all beta
    V_betas <- diag(c(V0, c(gammas * c2 * tau2 + (1 - gammas) * tau2)))
    D_beta <- solve(t(X) %*% X / sigsq + solve(V_betas))
    d_beta <- t(X) %*% y / sigsq
    H_beta <- t(chol(D_beta))
    beta <- D_beta %*% d_beta + H_beta %*% matrix(rnorm(kk + 1))
    # print(beta)
    # Conditional for all gamma
    for (j in 1:kk) {
        numerator <- p * dnorm(beta[j + 1], mean = 0, sd = sqrt(c2 * tau2))
        denominator <- numerator + (1 - p) * dnorm(beta[j + 1], mean = 0, sd = sqrt(tau2))
        prob <- numerator / denominator
        gammas[j, 1] <- ifelse(runif(1) < prob, 1, 0)
    }
    # print(gammas)
    # Conditional for sigma_sq
    e <- y - X %*% beta
    sigma_sq <- 1 / rgamma(1, shape = nn / 2 + a, scale = solve(1 / b + 0.5 * t(e) %*% e))
    # print(sigma_sq)

    if (i %% (niter / 10) == 0) {
        print(paste("loops... ", i))
    }

    if (i > burn) {
        nsim <- i - burn
        store_beta[, nsim] <- beta
        store_sigma[, nsim] <- sigma_sq
        store_gamma[, nsim] <- gammas
    }
}

apply(store_beta, MARGIN = 1, FUN = mean)
mean(store_sigma, na.rm = T)
apply(store_gamma, MARGIN = 1, FUN = mean)
