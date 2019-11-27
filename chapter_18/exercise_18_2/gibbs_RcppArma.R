#--------------------------------------------------------------------------------------------------
# R Code for Exercise 18.2
# Chan, J., Koop, G., Poirier, D.J. and Tobias, J.L. (2019).
# Bayesian Econometric Methods (2nd edition).
# Cambridge: Cambridge University Press.
#--------------------------------------------------------------------------------------------------

# run time: 8 sec when using the following compiler options specified in 
# ~/.R/Makevars options on a MacBookPro 10.14.6
# CC=/usr/local/bin/gcc-9
# CXX=/usr/local/bin/g++-9
# CXX1X=/usr/local/bin/g++-9
# CXX11=/usr/local/bin/g++-9
# SHLIB_CXXLD=/usr/local/bin/g++-9
# FC=/usr/local/bin/gfortran-9
# F77=/usr/local/bin/gfortran-9
# MAKE=make -j11
#  
# SHLIB_OPENMP_CFLAGS=-fopenmp
# SHLIB_OPENMP_CXXFLAGS=-fopenmp
# SHLIB_OPENMP_FCFLAGS=-fopenmp
# SHLIB_OPENMP_FFLAGS=-fopenmp
# PKG_CFLAGS=-fopenmp
# PKG_LIBS=-lgomp

library(Rcpp)
library(Matrix)
library(dplyr)
library(ggplot2)

rm(list = ls())

nsim <- 100000
nburn <- 2000
total_runs <- nsim + nburn

data <-
    read.csv("./chapter_18/exercise_18_2/usgdp.csv", header = FALSE)
y <- ts(
    data = data,
    start = c(1959, 1),
    freq = 4
)
y <- c(log(y)) * 100
y <- as.matrix(y)
colnames(y) <- NULL


sourceCpp("./chapter_18/exercise_18_2/gibbs_RcppArma.cpp",
          showOutput = TRUE,
          rebuild = TRUE,
          #dryRun = TRUE,
          verbose = TRUE
)

now <- Sys.time()
return_list <- gibbsC(nsim = nsim, nburn = nburn, y = y)
print(Sys.time() - now)

# Results #----------------------------------------------------------------------------------------
tau_hat <- rowMeans(return_list[[1]])
phi_hat <- rowMeans(return_list[[2]])
gamma_hat <- rowMeans(return_list[[3]])
sigma2_tau_hat <- mean(return_list[[4]])
sigma2_c_hat <- mean(return_list[[5]])
theta_hat <- matrix(c(c(phi_hat), sigma2_c_hat, sigma2_tau_hat, c(gamma_hat)))

cc <- ts(y - tau_hat, start = c(1949, 1), frequency = 4)

timetk::tk_tbl(data = cc, rename_index = "time") %>%
    rename(c = "Series 1") %>%
    ggplot(aes(x = time, y = c)) +
    geom_line()

print(theta_hat)
