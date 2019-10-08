# This function runs the Gibbs sampler 
# using the Rat growth data. 
#randn('seed',sum(100*clock));
rm(list=ls())
set.seed(123)
rats <- read.csv("./chapter_13/exercise_13_4/rats.txt", sep = "", header=FALSE)
rats <- as.matrix(rats)
dimnames(rats) <- NULL
#rats <- rats[,2:6]

nn <- nrow(rats)
kk <- ncol(rats)

xuse = cbind(1, c(8, 15, 22, 29, 36))

# -----------------------------
# Declare the prior values
# -----------------------------
a <- 3
b <- 1/(2*20)
    # this chooses the prior mean for 
    # sigma^2 equal to 20 with std
    # also equal to 20

eta <- matrix(c(100, 15), nrow=2)
C <- diag(2)
diag(C) <- c(40^2, 10^2)
        # this chooses the prior to center
        # the weight at birth at 100, and 
        # the growth rate at 15. C makes 
        # these choices reasonably diffuse

rho <- 5
R <-  diag(2)
diag(R) <- c(10^2, .5^2)
        #this prior specifies some degree 
        #of variation across rats, and does 
        #not restrict them to have equal 
        #birth weights and growth rates. 

        # for i = 1:1000;
        #     tempp = wish_rnd(inv(rho*R),rho);
        #     tempp2 = inv(tempp);
        #     std1(i,1) = sqrt(tempp2(1,1));
        #     std2(i,1) = sqrt(tempp2(2,2));
        # end;
        #  hist(std1,50);
        #  hist(std2,50);
        
iter <-  1e4 # declare the number of iterations
burn <-  5e2 # declare the length of the burn-in

#---------------------------------
# Set lengths of parmaeter vectors
#---------------------------------
sigma2 <-  matrix(0, nrow=iter, ncol=1)
theta0 <-  matrix(0, nrow=2,ncol=iter)
theta_int <-  matrix(0, nrow=30,ncol=iter)
theta_rate <-  matrix(0, nrow=30,ncol=iter)
invSigma <-  array(0, dim=c(2,2,iter))
var_int <-  matrix(0, nrow=iter, ncol=1)
var_rate <-  matrix(0, nrow=iter, ncol=1)
correl <-  matrix(0, nrow=iter, ncol=1)
        
#----------------------
# Set initial conditions
#----------------------
sigma2[1] <-  20
invSigma[,,1] <-  solve(matrix(c(100, 0, 0, 1), nrow=2,ncol=2))
theta0[,1] <-  matrix(c(100, 10), nrow=2)

# Begin the Gibbs Sampler
for (i in 2:iter) {
  
    total_resid <-  0
    # Do conditional for all theta_i
    for (j in 1:30) {

        Dtemp <- solve(t(xuse)%*%xuse/sigma2[i-1] + invSigma[,,i-1])
        dtemp <- t(xuse)%*%t(rats[j,2:kk, drop=FALSE])/sigma2[i-1] + invSigma[,,i-1]%*%theta0[,i-1]
        H <- chol(Dtemp);
        theta_temp = Dtemp%*%dtemp + t(H)%*%matrix(rnorm(2), nrow=2)
        theta_int[j,i] <- theta_temp[1,1]
        theta_rate[j,i] <- theta_temp[2,1]

        # use this later for sigma^2 conditional
        resids <- t(rats[j,2:kk, drop=FALSE]) - xuse%*%c(theta_int[j,i], theta_rate[j,i])
        tempp <-  sum(resids^2)
        total_resid <-  total_resid+tempp
    }
    thetabar = matrix(c(mean(theta_int[,i]), mean(theta_rate[,i])), nrow=2)

# Do conditional for theta0
Dtemp <-  solve(30*invSigma[,,i-1] + solve(C))
dtemp <-  30*invSigma[,,i-1]%*%thetabar + solve(C)%*%eta
H <-  chol(Dtemp)
theta0[,i] <- Dtemp%*%dtemp + t(H)%*%matrix(rnorm(2), nrow=2)

# Do Conditional for sigma2
sigma2[i,1] = 1/rgamma(n=1,shape=(150/2)+a, scale=(solve( .5*total_resid + solve(b)))^-1)

# Do Conditional for Sigma^-1. 
tempp3 <- 0
for (jj in 1:30) {
tempp1 <- matrix(c(theta_int[jj,i], theta_rate[jj,i]), nrow=2) - matrix(theta0[,i], nrow=2)
tempp2 <- tempp1%*%t(tempp1)
tempp3 <- tempp3 + tempp2
}

invSigma[,,i] <- rWishart(n=1, df = 30+rho, Sigma = solve(tempp3[,,drop=TRUE] + rho*R))

Sigma <- solve(invSigma[,,i])
var_int[i,1] <- Sigma[1,1]
var_rate[i,1] <- Sigma[2,2]
correl[i,1] <- Sigma[1,2]/(sqrt(Sigma[1,1])*sqrt(Sigma[2,2]))
}

int_use <- theta_int[,(burn+1):iter]
rate_use <- theta_rate[,(burn+1):iter]

mu_ints <-  mean(t(int_use))
mu_rates <-  mean(t(rate_use))


#save ratgibbs mu_ints mu_rates;

print('Posterior means of parameters')
print('Order: common intercept, common rate, variance_intercept, variance_rate, correlation')
print('intercept(10), intercept(25), rate(10), rate(25)')
parms <-  rbind(theta0[,(burn+1):iter],
var_int[(burn+1):iter],
var_rate[(burn+1):iter],
correl[(burn+1):iter],
theta_int[10,(burn+1):iter],
theta_int[25,(burn+1):iter],
theta_rate[10,(burn+1):iter],
theta_rate[25,(burn+1):iter]
)
resu <- data.frame('Post.mean'=rowMeans(parms),
                   'Post.std' = apply(parms, 1, sd),
                   '10th percentile' = apply(parms, 1, quantile, .10),
                    '90th percentile' = apply(parms, 1, quantile, .90))
rownames(resu) <- c("alpha_0", "beta_0", "sigma^2_alpha", "sigma^2_beta",
                    "rho_alpha_beta", "alpha_10", "alpha_25", "beta_10", "beta_25")
resu <- round(resu, 2)
print(resu)




