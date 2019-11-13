#--------------------------------------------------------------------------------------------------
# R Code for Exercise 18.2
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
#--------------------------------------------------------------------------------------------------

library(Rcpp)
rm(list = ls())

nsim <- 10000
nburn <- 1000
total_runs <- nsim+nburn

data <- read.csv("./chapter_18/exercise_18_2/usgdp.csv", header = FALSE)
y <- ts(data = data, start = c(1959, 1), freq = 4)
y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL

sourceCpp("./chapter_18/exercise_18_2/gibbs_RcppArma.cpp")
now <- Sys.time()
return_list <- gibbsC(nsim = nsim, nburn = nburn, y = y)
Sys.time()-now


# Results #----------------------------------------------------------------------------------------
tau_hat <- rowMeans(return_list[[1]])
phi_hat <- rowMeans(return_list[[2]])
gamma_hat <- rowMeans(return_list[[3]])
sigma2_tau_hat <- mean(return_list[[4]])
sigma2_c_hat <- mean(return_list[[5]])

cc <- ts(y-tau_hat, start = c(1949, 1), frequency = 4)



