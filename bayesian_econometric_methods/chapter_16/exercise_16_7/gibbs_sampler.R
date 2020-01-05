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
Sigma_sim <- matrix()

Sigma_sim <- diag(5)
Sigma_sim[lower.tri(Sigma_sim)] <- c(.4, .4, 0, .6, .7, 0, .3, 0, .3, 0)
Sigma_sim[upper.tri(Sigma_sim)] <- t(Sigma_sim)[upper.tri(Sigma_sim)]
L <- t(chol(Sigma_sim))
X_sim <- cbind(1, t(L %*% matrix(rnorm(5000), nrow = 5)))
beta_sim <- c(2, 25, -.4, -6, 0, 0)

y_sim <- X_sim %*% beta_sim + .2 * rnorm(1000)


nn <- nrow(y_sim)
kk <- 5

# xuse <- cbind(1, c(8, 15, 22, 29, 36))

niter <- 1e5 # declare the number of iterations
burn <- 2e2 # declare the length of the burn-in


# -----------------------------
# Declare the prior values
# -----------------------------
# prior for gamma_i
p <- 0.5
gammas <- matrix(1, nrow = 5)

# prior for beta
mu <- 0
tau2 <- 1e-8
c2 <- 9 / tau2
v0 <- 10^2
V_beta <- diag(kk + 1)
diag(V)[1] <- v0
diag(V)[2:(kk + 1)] <- gamma * matrix(rnorm(5, 0, tau^2)) + (1 - gamma) * matrix(rnorm(5, 0, (cc * tau)^2))


# prior for sigma^2
a <- 4
b <- 10

X <- X_sim
y <- y_sim
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
e <- (y - X %*% beta_hat)
ss <- t(e) %*% e
Sigma <- c(ss) * diag(kk + 1)



#---------------------------------
# Set lengths of parameter vectors
#---------------------------------
store_beta <- matrix(0, nrow = kk + 1, ncol = niter)
store_sigma <- matrix(0, nrow = 1, ncol = niter)
store_gamma <- matrix(0, nrow = kk, ncol = niter)

# sigma2 <- matrix(0, nrow = iter, ncol = 1)
# theta0 <- matrix(0, nrow = 2, ncol = iter)
# theta_int <- matrix(0, nrow = 30, ncol = iter)
# theta_rate <- matrix(0, nrow = 30, ncol = iter)
# invSigma <- array(0, dim = c(2, 2, iter))
# var_int <- matrix(0, nrow = iter, ncol = 1)
# var_rate <- matrix(0, nrow = iter, ncol = 1)
# correl <- matrix(0, nrow = iter, ncol = 1)

#----------------------
# Set initial conditions
#----------------------
Sigma <- c(ss) * diag(nn)
invSigma <- solve(Sigma)
beta <- matrix(0, nrow = kk + 1)
gamma <- ifelse(runif(5) > p, 1, 0)


#----------------------
# Begin the Gibbs Sampler
#----------------------
for (i in 1:(niter + burn)) {
    # Conditional for all beta
    XSig <- t(X) %*% invSigma
    Dtemp <- solve(XSig %*% X + solve(V_beta)) # / sigma2[i - 1] + invSigma[, , i - 1])
    dtemp <- (XSig %*% y) # t(xuse) %*% t(rats[j, 2:kk, drop = FALSE]) / sigma2[i - 1] + invSigma[, , i - 1] %*% theta0[, i - 1]
    H <- chol(Dtemp)
    beta <- Dtemp %*% dtemp + t(H) %*% matrix(rnorm(kk + 1))
    # print(beta)
    # Conditional for sigma_sq
    e <- y - X %*% beta
    # tempp <- sum(resids^2)
    # total_resid <- total_resid + tempp
    a_gamma <- nn / 2 + a
    b_gamma <- solve(1 / b + 0.5 * t(e) %*% e)
    sigma_sq <- 1 / rgamma(1, a_gamma, solve(b_gamma))
    diag(invSigma) <- sigma_sq^-1
    # print(sigma_sq)
    # print(invSigma)
    # Conditional for all gamma
    phi_p <- p / (cc * tau) * exp(-beta[-1, ]^2 / (2 * (cc * tau)^2))
    phi_1p <- (1 - p) / tau * exp(-beta[-1, ]^2 / (2 * (cc * tau)^2))
    prop_suc <- phi_p / (phi_p + phi_1p)
    for (i in 1:kk) {
        gamma <- ifelse(runif(5) > prop_suc, 1, 0)
    }
    # print(gamma)
    diag(V_beta)[2:(1 + kk)] <- c(gamma * tau^2 + (1 - gamma) * (tau * cc)^2)
    # V_beta <- diag()
    if (i %% 1000 == 0) {
        print(paste("loops... ", i))
    }

    if (niter > burn) {
          nsim <- niter - burn
      }
    store_beta[, nsim] <- beta
    store_sigma[, nsim] <- sigma_sq
    store_gamma[, nsim] <- gamma
}

round(apply(store_beta, MARGIN = 1, FUN = mean), 4)
round(mean(store_sigma, 4))
round(apply(store_gamma, MARGIN = 1, FUN = mean), 4)



int_use <- theta_int[, (burn + 1):iter]
rate_use <- theta_rate[, (burn + 1):iter]

mu_ints <- mean(t(int_use))
mu_rates <- mean(t(rate_use))


print("Posterior means of parameters")
print("Order: common intercept, common rate, variance_intercept, variance_rate, correlation")
print("intercept(10), intercept(25), rate(10), rate(25)")
parms <- rbind(
    theta0[, (burn + 1):iter],
    var_int[(burn + 1):iter],
    var_rate[(burn + 1):iter],
    correl[(burn + 1):iter],
    theta_int[10, (burn + 1):iter],
    theta_int[25, (burn + 1):iter],
    theta_rate[10, (burn + 1):iter],
    theta_rate[25, (burn + 1):iter]
)
resu <- data.frame(
    "Post.mean" = rowMeans(parms),
    "Post.std" = apply(parms, 1, sd),
    "10th percentile" = apply(parms, 1, quantile, .10),
    "90th percentile" = apply(parms, 1, quantile, .90)
)
rownames(resu) <- c(
    "alpha_0", "beta_0", "sigma^2_alpha", "sigma^2_beta",
    "rho_alpha_beta", "alpha_10", "alpha_25", "beta_10", "beta_25"
)
resu <- round(resu, 2)
print(resu)
