#------------------------------------------------------------------------------
# R Code for a fast Kalman Filter based on 
# Chan, J.C.C. and Jeliazkov, I. (2009). Efficient Simulation and
# Integrated Likelihood Estimation in State Space Models,
# International Journal of Mathematical Modelling and Numerical
# Optimisation, 1, 101-120.
#------------------------------------------------------------------------------

library(Matrix)
library(SuppDists)
library(RcppZiggurat)
library(Rcpp)


# TVP - VAR #---------------------------------------------------------------
# VAR:
# y[t] = mu[t] + Gamma[t] * y[t-1] + epsilon[t],   epsilon[t] ~ N(0, Omega_11)
# 
# State Space Form
# y[t] = X[t] * beta[t] + epsilon[t]   epsilon[t] ~ N(0, Omega_11)
# beta[t]=beta[t-1] + nu[t]     nu[t] ~ N(0, Omega_22)
#
# beta[t] = vec(mu[t], Gamma[t]')'
# 

rm(list = ls())
now <- Sys.time()

set.seed(123)

nsim <- 10
nburn <- 10
total_runs <- nsim + nburn

data <- read.csv("./papers/Chan_Jeliazkov_2009/USdata.csv", header = FALSE)
colnames(data) <- c("gpd_growth", "unemp", "inter", "inf")
y <- ts(data = data,
        start = c(1948, 1),
        freq = 4)

# ggplot2::ggplot(reshape2::melt(y), ggplot2::aes(x = Var1, y=value, col=Var2))+
#     ggplot2::geom_line()

#y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL

Y <- c(t(y))

TT <- nrow(y)
nn <- ncol(y)
qq <- nn*(nn+1)

# priors #---------------------------------------------------------------------



# Omega_11
nu_1 <- nn + 3 # nn?
S_1 <- diag(nn)

# Omega_22
DD <- 5 * diag(qq)
nu_2 <- rep(6, qq)
S_2 <- rep(0.01, qq)


# initial values #-------------------------------------------------------------

# Omega_11
Omega11 <- Matrix(cov(y))
Omega11_inv <- Matrix(solve(Omega11, diag(ncol(Omega11))))

# H
diags_H <- list(rep(1, TT*qq), rep(1, (TT - 1)*qq))
H <- bandSparse(n = TT*qq,
               k = c(0, -qq),
               diag = diags_H,
               symm = FALSE
               )

# S
DD_inv <- solve(DD)
Omega22 <- 0.01*diag(qq) # initial values for variance of betas?
Omega22_inv <- solve(Omega22)
S <- bdiag(replicate((TT-1), Omega22, simplify = F))
S <- bdiag(DD, S)  

S_inv <- bdiag(replicate((TT-1), Omega22_inv, simplify = F))
S_inv <-  bdiag(DD_inv, S_inv)  

# bigG = X
G <- bdiag(replicate(qq, matrix(c(1, y[1,]), ncol=nn+1), simplify = F))

G_list <- vector("list", TT)
for (i in 1:TT){
G_list[[i]] <- kronecker(diag(4), matrix(c(1, y[i,]), ncol=nn+1))
}
G <- bdiag(G_list)

# store
store_beta <- array(NA, dim=c(nsim, qq))
store_Omega11 <- array(NA, dim=c(nn, nn, nsim))
store_Omega22 <- array(NA, dim=c(nsim, qq))

# Gibbs Sampler -----------------------------------------------------------
for (ii in 1:total_runs) {

# beta --------------------------------------------------------------------
    S_inv <- bdiag(replicate((TT-1), Omega22_inv, simplify = F))
    S_inv <-  bdiag(DD_inv, S_inv)  
    K <- crossprod(H, S_inv) %*% H
    G_Omega11_inv <- crossprod(G, kronecker(diag(TT), Omega11_inv))
    G_Omega11_inv_G <- G_Omega11_inv %*% G
    P <- K + G_Omega11_inv_G
    L <- chol(P)
    
    beta_hat <- backsolve(L, forwardsolve(L, G_Omega11_inv  %*% Y, upper.tri = FALSE), upper.tri = TRUE,
              transpose = TRUE)
    beta <- beta_hat + solve(L, rnorm(n = length(beta_hat)))
    
# Omega_11 ----------------------------------------------------------------
    # v_1^0 & S_1^0 are parameters of the prior
    # y_t is observed
    # beta is sampled above
    # X_t needs to be designed
    e1 <- matrix(Y-G%*%beta, nrow=nn, ncol = TT)
    S_IW <- S_1 + tcrossprod(e1) 
    Omega11 <- 1/rWishart(n = 1, df = nu_1 + TT, Sigma = solve(S_IW, diag(nn)))[,,1]
    Omega11_inv <- solve(Omega11, diag(nn))

    
# Omega_22 |---------------------------------------------------------------
    Beta <- matrix(beta, ncol=qq)
    beta_D_sq <- colSums(diff(Beta)**2)
    for(i in seq_along( beta_D_sq)){
        
     diag(Omega22)[i] <- 1/rgamma(n=1,
                 shape = nu_2[i]+TT-1,
                 scale = 2/(S_2[i]+beta_D_sq[i])
                     )
    } 
    S_inv <- bdiag(replicate((TT-1), Omega22_inv, simplify = F))
    S_inv <-  bdiag(DD_inv, S_inv)  
    
    
# store |------------------------------------------------------------------
    if (ii > nburn) {
        run <- ii - nburn
        store_beta[run, ] <- beta
        store_Omega11[ , , run] <- Omega11
        store_Omega22[run, ] <- diag(Omega22)
    }
}
print(Sys.time() - now)


# Results #--------------------------------------------------------------------
tau_hat <- rowMeans(store_tau)
theta_hat <- rowMeans(store_theta)
print(theta_hat)

cc <- ts(y - tau_hat, start = c(1949, 1), frequency = 4)


timetk::tk_tbl(data = cc, rename_index = "time") %>%
    rename(c = "Series 1") %>%
    ggplot2::ggplot(aes(x = time, y = c)) +
    geom_line()
