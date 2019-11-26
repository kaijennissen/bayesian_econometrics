# TVP_VAR_DPS_DMA_sim.m - Forecasting with Large TVP-VAR using forgetting factors
# MULTIPLE MODEL CASE / DYNAMIC PRIOR SELECTION (DPS) AND DYNAMIC MODEL
# AVERAGING (DMA)
#-------------------------------------------------------------------------------
# The model is:
#
#	 y[t] = theta[t] x[t] + e[t]
#	 theta[t] = theta[t-1] + u[t]
#
# where x[t] = I x (y[t-1],...,y[t-p]) (Kronecker product), and e[t]~N(0,V[t])
# and u[t]~N(0,Q[t]).
#
# Additionally:
#
#  V[t] = kappa V[t-1] + (1-kappa) e[t-1]e[t-1]'
#  Q[t] = (1 - 1/lambda) S[t-1|t-1]
#
# This code estimates lambda and allows it to be time-varying. The specification is:
#
#  lambda[t] = lambda[min] + (1-lambda[min]) LL^(e[t]e[t]')
#
#-------------------------------------------------------------------------------
#  - This code allows to calculate ONLY iterated forecasts
#  - This code does predictive simulation by fixing theta[T+1],
#  theta[T+2],... from Normals with mean theta[T] (no drifting parameters
#  out-of-sample).
#  - This code does "online" forecasting, i.e. the Minnesota prior should not be
#  dependent on the data, so that the Kalman filter runs once for 1:T.
#-------------------------------------------------------------------------------

library(R.matlab)
rm(list = ls())

# Add path of data and functions
data_path <- "data"
func_path <- "functions"

#-------------------------------PRELIMINARIES--------------------------------------
# Choose grids for major tuning parameters
gamma  <-  c(1e-10, 1e-5, 0.001, 0.01, 0.05, 0.1, 1, 5)
lambda <- 0.99
kappa  <-  0.94

eta <-
    0.99   # Forgetting factor for DPS (dynamic prior selection) and DMA

# Please choose:
p <-  2            # p is number of lags in the VAR part
N <-  19            # Number of cross-sections (countries)

prior <- 1        # 1: Use Koop-type Minnesota prior
# 2: Use Litterman-type Minnesota prior

# Variables not to include in DMA
n_DMA <-  3      # Variables always included
varsN <-  0 + 1  # Additional variables in model averaging

# Forecasting
nfore <-  12       # Forecast horizon (note: forecasts are iterated)
t0 <-
    "2005M12"    # Set last observation of initial estimation period
nsim <-
    500      # Number of times to simulate from the predictive density

# Choose which results to print
# NOTE: CHOOSE ONLY 0/1 (FOR NO/YES) VALUES!
print_fore <- 1           # summary of forecasting results
print_coefficients <-
    0   # plot volatilities and lambda_t (but not theta_t which is huge)
print_pred <- 0           # plot predictive likelihoods over time
print_Min <- 0            # print the Minnesota prior over time
#----------------------------------LOAD DATA----------------------------------------
# ========| LOAD DATA
readMat("./bvars/LARGE_PVAR/data/PVARdataEA.mat")
Y <- Ynew
Y <- Y[2:end,]
yearlab <- yearlab[(p + 2):end]

T_thres <-
    find(strcmp(yearlab, t0) == 1)  # Convert t0 to numeric value

for (country_index in 1:N) {
    Y1 <- cell(varsN, 1)
    M <- matrix(0, nrows = varsN, ncols = 1)
    
    for (ss in 1:varsN) {
        G <-
            ss + n_DMA    # Number of macro fundamentals + the exchange rate
        SS1 <- (1:G)
        SS <- vector(mode = "numeric")
        GG <- SS1 + 7 * (country_index - 1)
        SS <- [SS, GG]
        Y1[ss, 1] <- Y[:, SS]
        M[ss, 1] <- length(SS) # M is the dimensionality of Y
    }
    t <- size(Y1[1, 1], 1)
    
    nos <- varsN
    
    # Inflation is the variable of interest for forecasting
    nfocus = 1
    
    # ===================================| VAR EQUATION |==============================
    # Generate lagged Y matrix. This will be part of the X matrix
    x_t = cell(nos, 1)
    x_f = cell(nos, 1)
    y_t = cell(nos, 1)
    K = matrix(0, nrow=nos, ncol=1)
    
    for (ss in 1:nos) {
        ylag = mlag2(Y1[ss, 1], p)
        ylag = ylag(p + 1:end,:)
        [temp, kk] = create_RHS(ylag, M(ss), p, t)
        x_t[ss, 1] = [ones(size(ylag, 1), 1) ylag]
        K(ss, 1) = kk
        x_f[ss, 1] = ylag
        y_t[ss, 1] = Y1[ss, 1](p + 1:end,:)
    }
    
    yearlab = yearlab(p + 1:end)
    # Time series observations
    t = size(y_t[1, 1], 1) ##ok<*NASGU>
    
    #----------------------------PRELIMINARIES---------------------------------
    
    anumber = t - T_thres + 1
    y_fore = cell(nos, 1)
    
    if (country_index == 1) {
        LOG_PL_VAR = matrix(0, nrow = anumber, ncol = nfore)
        MSFE_VAR = zeros(anumber, N * nfocus, nfore)
        MAFE_VAR = zeros(anumber, N * nfocus, nfore)
        logpl_VAR = zeros(anumber, N * nfocus, nfore)
    }
    
    offset = 1e-9  # just a constant for numerical stability
    
    #----------------------------- END OF PRELIMINARIES ---------------------------
    
    #======================= BEGIN KALMAN FILTER ESTIMATION =======================
    
    for (irep in T_thres:t) {
        if (mod(irep, ceil(t. / 40)) == 0) {
            disp t([num2str(100 * (irep / t))]) # completed'
            toc
        }
        
        beta_OLS = (
            x_t[1, 1](1:irep,:)'*x_t[1, 1](1:irep,:))\(x_t[1, 1](1:irep,:)' * y_t{
                1, 1
            }(1:irep,:)
        )
        sigma_OLS = (y_t[1, 1]}(1:irep,:) - x_t[1, 1](1:irep,:) * beta_OLS)'*(y_t[1, 1](1:irep,:) - x_t[1, 1](1:irep,:)*beta_OLS)/(irep-M(1))

        if (irep >= T_thres){

            # Start predictive simulation
            chol_S = chol(sigma_OLS)

            for (sim in 1:nsim){

                Y_hat = 0
                # Now create forecast for h=1
                X_FORE = [1 y_t{ss,1}(irep,:) x_f{ss,1}(irep,1:M(ss)*(p-1))]
                Y_hat = X_FORE*beta_OLS + randn(1,M(ss))*chol_S
                y_fore{ss,1}(1,:,sim) = Y_hat
                # Now do forecasts for h>1

                for (ii in 1:nfore-1) {
                    if (ii <= p) {   # if h<=p (number of lags)
                        X_new_temp = [1 Y_hat X_FORE(:,2:M(ss)*(p-ii)+1)]
                        Y_temp = X_new_temp*beta_OLS + randn(1,M(ss))*chol_S
                        Y_hat = [Y_temp Y_hat]
                    } else { # if h>p (number of lags)
                        X_new_temp = [1 Y_hat(:,1:M(ss)*p)]
                        Y_temp = X_new_temp*beta_OLS + randn(1,M(ss))*chol_S
                        Y_hat = [Y_temp Y_hat]
                }

                    # This cell array saves the draws from all forecast
                    # horizons, from all model sizes.
                    y_fore{ss,1}(ii+1,:,sim) = Y_temp

                }
            } # End predictive simulation

            # Find "observed" out-of-sample data for MSFE and MAFE calculations
            if (irep <= t-nfore) {
                Yraw_f[1, 1] = y_t[1, 1](irep+1:irep+nfore,:,:) #Pseudo out-of-sample observations
            } else {
                Yraw_f[1, 1] =[y_t[1, 1](irep+1:t,:)  NaN(nfore-(t-irep),M(1))]
            }

            # Now we have the predictions for each model & the associated model
            # probabilities

            for (ii in 1:nfore){

                focus_vars = 1
                y_t_VAR(ii,:,irep-T_thres+1) = mean(y_fore{ss,1}(ii,focus_vars,:),3)
                variance_VAR = cov(squeeze(y_fore{ss,1}(ii,focus_vars,:))')

LOG_PL_VAR(irep-T_thres+1,ii) = log(mvnpdfs(Yraw_f[1, 1](ii,focus_vars)',y_t_VAR(ii,:,irep-T_thres+1)',variance_VAR) + offset)
MAFE_VAR(irep-T_thres+1,(country_index-1)*nfocus+1:country_index*nfocus,ii) = abs(Yraw_f[1, 1](ii,focus_vars) - squeeze(y_t_VAR(ii,:,irep-T_thres+1)))
MSFE_VAR(irep-T_thres+1,(country_index-1)*nfocus+1:country_index*nfocus,ii) = (Yraw_f[1, 1](ii,focus_vars) - squeeze(y_t_VAR(ii,:,irep-T_thres+1))).^2

MAFE_RW(irep-T_thres+1,(country_index-1)*nfocus+1:country_index*nfocus,ii) = abs(y_t[1, 1](irep,focus_vars)' - Yraw_f[1, 1](ii,focus_vars))
                MSFE_RW(irep-T_thres+1,(country_index-1)*nfocus+1:country_index*nfocus,ii) = (y_t[1, 1](irep,focus_vars)' - Yraw_f[1, 1](ii,focus_vars)).^2

j_in = 0

for (j in focus_vars){
    j_in = j_in + 1
    logpl_VAR(irep-T_thres+1,(country_index-1)*nfocus+1:country_index*nfocus,ii) = log(mvnpdfs(Yraw_f{ss}(ii,j)',y_t_VAR(ii,j_in,irep-T_thres+1)',variance_VAR(j_in,j_in)) + offset)
}
    }
}

}
#======================== END KALMAN FILTER ESTIMATION ========================

}

#===================| PRINT RESULTS |=========================

# format short g
# summary of forecasting results
if (print_fore == 1) {
    disp('MSFE for the key variables of interest')
    msfe_focus = mean(MSFE_VAR(1:end - 1,:, 1))
    for (iii = 2:nfore) {
        msfe_focus = [msfe_focus mean(MSFE_VAR(1:end - iii,:, iii))]
    }
    
    disp(['    Horizon    GDP       INFL      INTR']) ##ok<*NBRAK>
    disp([(1:nfore)' msfe_focus])

    disp('MAFE for the key variables of interest')
    mafe_focus = mean(MAFE_VAR(1:end-1,:,1))

    for (iii in 2:nfore) {
        mafe_focus=[mafe_focus mean(MAFE_VAR(1:end-iii,:,iii))]
    }

    disp(['    Horizon    GDP       INFL      INTR'])
    disp([(1:nfore)' mafe_focus])
    
    disp('sum of log pred likes for the key variables of interest individually')
    lpls_focus = sum(logpl_VAR(1:end - 1,:, 1))
    
    for (iii in 2:nfore) {
        lpls_focus = [lpls_focus sum(logpl_VAR(1:end - iii,:, iii))]
    }
    
    disp(['    Horizon    GDP      INFL      INTR'])
    disp([(1:nfore)' lpls_focus])
    disp('sum of log pred likes for the key variables of interest jointly')
    lpl_focus = sum(LOG_PL_VAR(1:end-1,1))

    for (iii in 2:nfore){
        lpl_focus=[lpl_focus sum(LOG_PL_VAR(1:end-iii,iii))]
    }

    disp(['    Horizon   TOTAL'])
    disp([(1:nfore)' lpl_focus])
    disp('                      ')
    disp('                      ')
    
}


# plot volatilities and lambda_t (but not theta_t which is huge)
if (print_coefficients == 1) {
    prob_pl = omega_update{
        1, 1
    }(index_DMA(:, 1))
    for (ss in 2:nos) {
        prob_pl = [prob_pl omega_update{
            ss, 1
        }(index_DMA(:, ss))] ##ok<*AGROW>
    }
    
    final_prob = 0 * prob_pl
    
    for (sss in 1:nos) {
        final_prob(:, sss) = prob_pl(:, sss). / sum(prob_pl, 2)
    }
    figure
    plot(yearlab, final_prob)
    warning off
    legend('small VAR', 'medium VAR', 'large VAR')
    title('Time-varying probabilities of small/medium/large VARs')
}

# plot predictive likelihoods over time
if (print_pred == 1) {
    figure
    subplot(2, 2, 1)
    plot(yearlab(T_thres + 1:end), cumsum(LOG_PL_DMS(1:end - 1, 1)))
    title('Cumulative sum of total predictive likelihood h=1')
    subplot(2, 2, 2)
    plot(yearlab(T_thres + 2:end), cumsum(LOG_PL_DMS(1:end - 2, 2)))
    title('Cumulative sum of total predictive likelihood h=2')
    subplot(2, 2, 3)
    plot(yearlab(T_thres + 3:end), cumsum(LOG_PL_DMS(1:end - 3, 3)))
    title('Cumulative sum of total predictive likelihood h=3')
    subplot(2, 2, 4)
    plot(yearlab(T_thres + 4:end), cumsum(LOG_PL_DMS(1:end - 4, 4)))
    title('Cumulative sum of total predictive likelihood h=4')
}
# print the Minnesota prior over time
if (print_Min == 1) {
    print(c(
        "=========MINNESOTA PRIOR==========",
        "\n",
        "Best gammas for each time period"
    ))
    vars_Min = t(T_thres:t)
    for (ss in 1:nos) {
        vars_Min = cbind(vars_Min  Minn_gamms[ss, 1])
    }
    lll = c("period", "smallVAR", "mediumVAR", "largeVAR")
    print(lll(1:nos + 1))
    print(vars_Min)
}

#save(sprintf('#s.mat','competing_forecasts'),'-mat')