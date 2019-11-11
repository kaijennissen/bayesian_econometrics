#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Ziggurat.h>
// [[Rcpp::depends(RcppArmadillo, Rcpp, RcppZiggurat)]]

using namespace Rcpp;
using namespace arma;
static Ziggurat::Ziggurat::Ziggurat zigg;

//[[Rcpp::export]]
Rcpp::NumericVector zrnorm(int n)
{
    Rcpp::NumericVector x(n);
    for (int i = 0; i < n; i++)
    {
        x[i] = zigg.norm();
    }
    return x;
}

// [[Rcpp::export]]
List gibbsC(int nsim, int nburn, mat y)
{

    List to_return(5);
    int TT = y.n_rows;

    // ----------------------------------------------------------------
    // Priors
    // ----------------------------------------------------------------
    // phi
    mat phi_0(2, 1);
    phi_0(0, 0) = 1.3;
    phi_0(1, 0) = -0.7;
    mat V_phi = eye(2, 2);
    mat invV_phi = inv(V_phi);

    // gamma
    mat gamma_0(2, 1);
    gamma_0.fill(750);
    mat V_gamma = 100 * eye(2, 2);
    mat invV_gamma = inv(V_gamma);

    // sigma2_tau
    double b_sigma2_tau = 0.01;

    // sigma2_c
    double ny_sigma2_c = 3;
    double S_sigma2_c = 2;

    // ----------------------------------------------------------------
    // Initial Values
    // ----------------------------------------------------------------
    mat phi(2, 1);
    phi(0, 0) = -1.34;
    phi(1, 0) = 0.7;
    double sigma2_c = .5;
    double sigma2_tau = .001;
    mat tau_0(2, 1);
    tau_0.fill(y(0));

    // H_2
    sp_mat H_2(TT, TT);
    for (int j = 0; j < TT; j++)
    {
        H_2(j, j) = 1;
        if (j > 0)
        {
            H_2(j, j - 1) = -2;
        }
        if (j > 1)
        {
            H_2(j, j - 2) = 1;
        }
    }
    mat invH_2 = inv(mat(H_2));
    sp_mat HH_2 = H_2.t() * H_2;

    // H_phi
    sp_mat H_phi(TT, TT);
    for (int j = 0; j < TT; j++)
    {
        H_phi(j, j) = 1;
        if (j > 0)
        {
            H_phi(j, j - 1) = phi(0);
        }
        if (j > 1)
        {
            H_phi(j, j - 2) = phi(1);
        }
    }
    sp_mat HH_phi = H_phi.t() * H_phi;

    // X_gamma
    mat X_gamma(TT, 2);
    for (int j = 0; j < TT; j++)
    {
        X_gamma(j, 0) = j + 2;
        X_gamma(j, 1) = -(j + 1);
    }

    // ----------------------------------------------------------------
    // set up intermediate data structures and storage variables
    // ----------------------------------------------------------------

    // intermediate vars tau
    vec S_sigma2_c_temp(1, 1, fill::zeros);
    mat alpha_tau_tilde(TT, 1);
    mat alpha_tau(TT, 1);
    sp_mat K_tau(TT, TT);
    mat XX_tau(TT, 1);
    mat tau_hat(TT, 1);
    mat tau(TT, 1);

    // intermediate vars phi
    mat K_phi(size(invV_phi));
    mat XX(size(K_phi));
    mat phi_hat(size(phi_0));
    mat phic(size(phi_hat));
    mat c(size(y));
    mat X_phi(TT, 2, fill::zeros);

    // sigma2_tau
    int n_grid = 500;
    mat del_tau(TT + 1, 1);
    mat sigma2_tau_grid(n_grid, 1);
    mat lp_sigtau2(size(sigma2_tau_grid));
    mat p_sigtau2(size(lp_sigtau2));
    mat cdf_sigtau2(size(p_sigtau2));
    uvec idx;

    // gamma
    mat K_gamma(size(invV_gamma));
    mat XX_temp_gamma(size(K_gamma));
    mat gamma_hat(size(gamma_0));
    mat gamma(size(gamma_hat));

    // Store
    mat tau_store(y.n_rows, nsim, fill::zeros);
    mat phi_store(2, nsim, fill::zeros);
    mat sigma2_c_store(1, nsim, fill::zeros);
    mat sigma2_tau_store(1, nsim, fill::zeros);
    mat gamma_store(2, nsim, fill::zeros);
    vec gamma_draw(1, 1, fill::zeros);

    for (int i = 0; i < nsim + nburn; i++)
    {

        // ------------------------------------------------------------
        // draw from conditional for tau
        // ------------------------------------------------------------
        alpha_tau_tilde = join_cols(vec{2 * tau_0(1) - tau_0(0), -tau_0(1)}, zeros(TT - 2));
        alpha_tau = mat(H_2).i() * alpha_tau_tilde;
        K_tau = HH_phi.t() / sigma2_c + HH_2 / sigma2_tau;
        XX_tau = HH_phi * y / sigma2_c + HH_2 * alpha_tau / sigma2_tau;
        tau_hat = solve(mat(K_tau), XX_tau);
        tau = tau_hat + inv(chol(mat(K_tau), "upper")) * randn(H_2.n_cols, 1);

        // ------------------------------------------------------------
        // draw from conditional for phi
        // ------------------------------------------------------------
        c = y - tau;

        for (unsigned int j = 1; j < TT; j++)
        {
            X_phi(j, 0) = c(j - 1, 0);
        }
        for (unsigned int j = 2; j < TT; j++)
        {
            X_phi(j, 1) = c(j - 2, 0);
        }
        K_phi = invV_phi + X_phi.t() * X_phi / sigma2_c;
        XX = invV_phi * phi_0 + X_phi.t() * c / sigma2_c;
        phi_hat = solve(K_phi, XX);
        phic = phi_hat + inv(chol(K_phi)) * randn(2, 1);
        if (phic(0) + phic(1) < .99 and phic(1) - phic(0) < .99 and phic(1) > -.99)
        {
            phi = phic;
            // calculate H_phi in every iteration
            for (int j = 0; j < TT; j++)
            {
                H_phi(j, j) = 1;
                if (j > 0)
                {
                    H_phi(j, j - 1) = phi(0);
                }
                if (j > 1)
                {
                    H_phi(j, j - 2) = phi(1);
                }
            }
            HH_phi = H_phi.t() * H_phi;
        }

        // ------------------------------------------------------------
        // draw from condition for sigma^2_c
        // ------------------------------------------------------------
        S_sigma2_c_temp = S_sigma2_c + 0.5 * (y - tau).t() * (H_phi * (y - tau));
        sigma2_c = 1 / arma::randg<double>(distr_param(ny_sigma2_c + TT / 2, 1 / as_scalar(S_sigma2_c_temp)));

        // ------------------------------------------------------------
        // draw from conditional for sigma^2_tau
        // ------------------------------------------------------------
        for (int j = 0; j < TT + 1; j++)
        {
            if (j > 1)
            {
                del_tau(j, 0) = tau(j - 1) - tau(j - 2);
            }
            else if (j == 0)
            {
                del_tau(j, 0) = tau_0(1) - tau_0(0);
            }
            else if (j == 1)
            {
                del_tau(j, 0) = tau(0) - tau_0(1);
            }
        }
        sigma2_tau_grid = linspace(as_scalar(randu(1)) / 100, b_sigma2_tau - as_scalar(randu(1)) / 100, n_grid);
        for (int j = 0; j < n_grid; j++)
        {
            lp_sigtau2(j, 0) = as_scalar(-log(sigma2_tau_grid(j)) * TT / 2 - sum(square(diff(del_tau))) / (2 * sigma2_tau_grid(j)));
        }
        p_sigtau2 = arma::exp(lp_sigtau2 - as_scalar(max(lp_sigtau2))) / as_scalar(sum(arma::exp(lp_sigtau2 - as_scalar(max(lp_sigtau2)))));
        cdf_sigtau2 = cumsum(p_sigtau2);
        idx = find(randu<double>() < cdf_sigtau2, 1);
        sigma2_tau = as_scalar(sigma2_tau_grid(idx));

        // ------------------------------------------------------------
        // draw from conditional for gamma
        // ------------------------------------------------------------
        mat K_gamma = invV_gamma + X_gamma.t() * (HH_2 * X_gamma) / as_scalar(sigma2_tau);
        mat XX_temp_gamma = invV_gamma * gamma_0 + X_gamma.t() * (HH_2 * tau) / as_scalar(sigma2_tau);
        mat gamma_hat = solve(K_gamma, XX_temp_gamma);
        mat gamma = gamma_hat + inv(chol(K_gamma)) * randn(2, 1);

        // ------------------------------------------------------------
        // Store draws after burn-in period
        // ------------------------------------------------------------
        if (i >= nburn)
        {
            tau_store.col(i - nburn) = tau;
            phi_store.col(i - nburn) = phi;
            sigma2_c_store.col(i - nburn) = sigma2_c;
            sigma2_tau_store.col(i - nburn) = sigma2_tau;
            gamma_store.col(i - nburn) = gamma;
        }
    }

    to_return(0) = tau_store;
    to_return(1) = phi_store;
    to_return(2) = sigma2_c_store;
    to_return(3) = sigma2_c_store;
    to_return(4) = gamma_store;

    return to_return;
}
