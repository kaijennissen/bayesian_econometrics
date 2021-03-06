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

#/usr/local/Cellar/r/3.6.1_1/lib/R/bin/R CMD SHLIB -o 'sourceCpp_14.so' --preclean  'gibbs_RcppArma.cpp'
#/usr/local/bin/g++-9 -I"/usr/local/Cellar/r/3.6.1_1/lib/R/include" -DNDEBUG -I../inst/include   -I"/usr/local/lib/R/3.6/site-library/Rcpp/include" -I"/usr/local/lib/R/3.6/site-library/RcppArmadillo/include" -I"/Users/kaijennissen/Intern/05_Projekte/bayesian_econometrics/chapter_18/exercise_18_2" -I/usr/local/opt/gettext/include -I/usr/local/opt/readline/include -I/usr/local/include  -fPIC  -g -O2  -c gibbs_RcppArma.cpp -o gibbs_RcppArma.o
#/usr/local/bin/g++-9 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/opt/gettext/lib -L/usr/local/opt/readline/lib -L/usr/local/lib -L/usr/local/Cellar/r/3.6.1_1/lib/R/lib -L/usr/local/opt/gettext/lib -L/usr/local/opt/readline/lib -L/usr/local/lib -o sourceCpp_14.so gibbs_RcppArma.o -lgomp -L/usr/local/Cellar/r/3.6.1_1/lib/R/lib -lR -lintl -Wl,-framework -Wl,CoreFoundation

#/usr/local/Cellar/r/3.6.1_1/lib/R/bin/R CMD SHLIB -o 'sourceCpp_15.so' --preclean  'gibbs_RcppArma.cpp'
#clang++ -std=gnu++11 -I"/usr/local/Cellar/r/3.6.1_1/lib/R/include" -DNDEBUG -I../inst/include   -I"/usr/local/lib/R/3.6/site-library/Rcpp/include" -I"/usr/local/lib/R/3.6/site-library/RcppArmadillo/include" -I"/Users/kaijennissen/Intern/05_Projekte/bayesian_econometrics/chapter_18/exercise_18_2" -I/usr/local/opt/gettext/include -I/usr/local/opt/readline/include -I/usr/local/include  -fPIC  -g -O2  -c gibbs_RcppArma.cpp -o gibbs_RcppArma.o
#clang++ -std=gnu++11 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/opt/gettext/lib -L/usr/local/opt/readline/lib -L/usr/local/lib -L/usr/local/Cellar/r/3.6.1_1/lib/R/lib -L/usr/local/opt/gettext/lib -L/usr/local/opt/readline/lib -L/usr/local/lib -o sourceCpp_15.so gibbs_RcppArma.o -L/usr/local/opt/openblas/lib -lopenblas -L/usr/local/opt/gcc/lib/gcc/9/gcc/x86_64-apple-darwin18/9.2.0 -L/usr/local/opt/gcc/lib/gcc/9 -lgfortran -lquadmath -lm -lgomp -L/usr/local/Cellar/r/3.6.1_1/lib/R/lib -lR -lintl -Wl,-framework -Wl,CoreFoundation


library(Rcpp)
library(Matrix)
library(dplyr)
library(ggplot2)

rm(list = ls())
create_custom_makevars <- function() {
    if (file.exists("~/.R/Makevars")) {
        file.rename(from = "~/.R/Makevars", to = "~/.R/Makevars_copy")
    }
    
    cat(" CC=/usr/local/bin/gcc-9" ,
        file = "~/.R/Makevars",
        sep = "\n")
    
    makevars_lines <- c(
        "CC=/usr/local/bin/gcc-9",
        "CC=/usr/local/bin/gcc-9",
        "CXX=/usr/local/bin/g++-9",
        "CXX1X=/usr/local/bin/g++-9",
        "CXX11=/usr/local/bin/g++-9",
        "SHLIB_CXXLD=/usr/local/bin/g++-9",
        "FC=/usr/local/bin/gfortran-9",
        "F77=/usr/local/bin/gfortran-9",
        "MAKE=make -j11",
        "SHLIB_OPENMP_CFLAGS=-fopenmp",
        "SHLIB_OPENMP_CXXFLAGS=-fopenmp",
        "SHLIB_OPENMP_FCFLAGS=-fopenmp",
        "SHLIB_OPENMP_FFLAGS=-fopenmp",
        #"PKG_CFLAGS=-fopenmp",
        "CXXFLAGS = -fopenmp -g -O2",
        #"PKG_LIBS=-lgomp"
        "PKG_LIBS=-lgomp"
    )
    
    for (makevars_line in makevars_lines) {
        cat(
            makevars_line,
            file = "~/.R/Makevars",
            sep = "\n",
            append = TRUE
        )
    }
    
    
    
    
}

remove_custom_makevars <- function() {
    if (file.exists("~/.R/Makevars_copy")) {
        if (file.exists("~/.R/Makevars")) {
            file.remove("~/.R/Makevars")
        }
        file.rename(from = "~/.R/Makevars_copy", to = "~/.R/Makevars")
    }
}


create_custom_makevars()
sourceCpp(
    "./bayesian_econometric_methods/chapter_18/exercise_18_2/gibbs_RcppArma.cpp",
    showOutput = TRUE,
    rebuild = TRUE,
    #dryRun = TRUE,
    verbose = TRUE
)
remove_custom_makevars()

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



now <- Sys.time()
return_list <- gibbsC(nsim = nsim, nburn = nburn, y = y)
print(Sys.time() - now)

# Results #----------------------------------------------------------------------------------------
tau_hat <- rowMeans(return_list[[1]])
phi_hat <- rowMeans(return_list[[2]])
gamma_hat <- rowMeans(return_list[[3]])
sigma2_tau_hat <- mean(return_list[[4]])
sigma2_c_hat <- mean(return_list[[5]])
theta_hat <-
    matrix(c(c(phi_hat), sigma2_c_hat, sigma2_tau_hat, c(gamma_hat)))

cc <- ts(y - tau_hat, start = c(1949, 1), frequency = 4)

timetk::tk_tbl(data = cc, rename_index = "time") %>%
    rename(c = "Series 1") %>%
    ggplot(aes(x = time, y = c)) +
    geom_line()

print(theta_hat)

