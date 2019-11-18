#include <iostream>
#include <armadillo>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>

using namespace arma;
using namespace std;

int main()
{

    string fileName = "usgdp.csv";
    int dPrec = 20;

    ifstream inputData;
    inputData.open(fileName);
    cout.precision(dPrec);
    mat outputMatrix;

    if (!inputData)
        return -1;
    string fileline, filecell;

    unsigned int prevNoOfCols = 0, noOfRows = 0, noOfCols = 0;

    while (getline(inputData, fileline))
    {
        noOfCols = 0;
        stringstream linestream(fileline);
        while (getline(linestream, filecell, ','))
        {
            try
            {
                stod(filecell);
            }
            catch (...)
            {
                return -1;
            }
            noOfCols++;
        }
        if (noOfRows++ == 0)
            prevNoOfCols = noOfCols;
        if (prevNoOfCols != noOfCols)
            return -1;
    }
    inputData.close();
    outputMatrix.resize(noOfRows, noOfCols);

    inputData.open(fileName);
    noOfRows = 0;
    while (getline(inputData, fileline))
    {
        noOfCols = 0;
        stringstream linestream(fileline);
        while (getline(linestream, filecell, ','))
        {
            outputMatrix(noOfRows, noOfCols++) = stod(filecell);
        }
        noOfRows++;
    }

    mat y = outputMatrix;
    y = log(y) * 100;

    int TT = y.n_rows;
    int nsim = 100000;
    int nburn = 20000;

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
    phi(0, 0) = 1.34;
    phi(1, 0) = -0.7;
    double sigma2_c = .5;
    double sigma2_tau = .001;
    mat gamma(2, 1);
    gamma.fill(y(0));

    // H_2
    umat H2_location(2, 3 * (TT - 1));
    vec H2_values(3 * (TT - 1));
    int i = 0;
    for (int j = 0; j < TT; j++)
    {
        H2_location(0, i) = j;
        H2_location(1, i) = j;
        H2_values(i) = 1;
        i++;
        if (j < TT - 1)
        {
            H2_location(0, i) = j + 1;
            H2_location(1, i) = j;
            H2_values(i) = -2;
            i++;
        }
        if (j < TT - 2)
        {
            H2_location(0, i) = j + 2;
            H2_location(1, i) = j;
            H2_values(i) = 1;
            i++;
        }
    }
    sp_mat H_2 = sp_mat(H2_location, H2_values);
    mat invH_2 = inv(mat(H_2));
    sp_mat HH_2 = H_2.t() * H_2;

    // H_phi
    uvec idx1(TT);
    uvec idx2(TT - 1);
    vec vals_phi_2(TT - 1);
    uvec idx3(TT - 2);
    vec vals_phi_3(TT - 2);
    umat Hphi_location(2, 3 * (TT - 1));
    vec Hphi_values(3 * (TT - 1));
    i = 0;
    for (int j = 0; j < TT; j++)
    {
        Hphi_location(0, i) = j;
        Hphi_location(1, i) = j;
        Hphi_values(i) = 1;
        idx1(j) = i;
        i++;
        if (j < TT - 1)
        {
            Hphi_location(0, i) = j + 1;
            Hphi_location(1, i) = j;
            Hphi_values(i) = -phi(0);
            idx2(j) = i;
            i++;
        }
        if (j < TT - 2)
        {
            Hphi_location(0, i) = j + 2;
            Hphi_location(1, i) = j;
            Hphi_values(i) = -phi(1);
            idx3(j) = i;
            i++;
        }
    }
    vals_phi_2.fill(-phi(0));
    vals_phi_3.fill(-phi(1));
    Hphi_values.elem(idx2) = vals_phi_2;
    Hphi_values.elem(idx3) = vals_phi_3;
    sp_mat H_phi = sp_mat(Hphi_location, Hphi_values);
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
    vec alpha_tau_tilde = zeros(TT);
    vec alpha_tau(TT);
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
    vec del_tau(TT + 1);
    vec sigma2_tau_grid(n_grid, 1);
    vec lp_sigtau2(size(sigma2_tau_grid));
    vec p_sigtau2(size(lp_sigtau2));
    vec cdf_sigtau2(size(p_sigtau2));
    uvec idx;

    // sigma2_c
    double ny_hat_sigma2_c = ny_sigma2_c + TT / 2;
    double S_hat_sigma2_c = 99.0;

    // gamma
    mat K_gamma(size(invV_gamma));
    mat XX_temp_gamma(size(K_gamma));
    mat gamma_hat(size(gamma_0));

    // // Store
    mat tau_store(y.n_rows, nsim, fill::zeros);
    mat phi_store(2, nsim, fill::zeros);
    vec sigma2_c_store(nsim, fill::zeros);
    vec sigma2_tau_store(nsim, fill::zeros);
    mat gamma_store(2, nsim, fill::zeros);
    vec gamma_draw(1, 1, fill::zeros);

    for (int i = 0; i < nsim + nburn; i++)
    {

        // ------------------------------------------------------------
        // draw from conditional for tau
        // ------------------------------------------------------------
        // gamma = (tau(0), tau(-1))'
        alpha_tau_tilde(0) = 2 * gamma(0) - gamma(1);
        alpha_tau_tilde(1) = -gamma(0);
        alpha_tau = invH_2 * alpha_tau_tilde;
        K_tau = HH_phi / sigma2_c + HH_2 / sigma2_tau;
        XX_tau = HH_phi * y / sigma2_c + HH_2 * alpha_tau / sigma2_tau;
        tau_hat = solve(mat(K_tau), XX_tau);
        tau = tau_hat + solve(chol(mat(K_tau), "upper"), eye(TT, TT)) * randn(H_2.n_cols, 1);

        // ------------------------------------------------------------
        // draw from conditional for phi
        // ------------------------------------------------------------
        c = y - tau;
        for (unsigned int j = 0; j < TT; j++)
        {
            if (j > 1)
            {
                X_phi(j, 0) = c(j - 1, 0);
                X_phi(j, 1) = c(j - 2, 0);
            }
            else if (j == 1)
            {
                X_phi(j, 0) = c(j - 1, 0);
                X_phi(j, 1) = 0;
            }
            else if (j == 0)
            {
                X_phi(j, 0) = 0;
                X_phi(j, 1) = 0;
            }
        }
        K_phi = invV_phi + X_phi.t() * X_phi / sigma2_c;
        XX = invV_phi * phi_0 + X_phi.t() * c / sigma2_c;
        phi_hat = solve(K_phi, XX);
        phic = phi_hat + inv(chol(K_phi, "upper")) * randn(2, 1);
        if (phic(0) + phic(1) < .99 and phic(1) - phic(0) < .99 and phic(1) > -.99)
        {
            phi = phic;
            // calculate H_phi in every iteration
            vals_phi_2.fill(-phi(0));
            vals_phi_3.fill(-phi(1));
            Hphi_values.elem(idx2) = vals_phi_2;
            Hphi_values.elem(idx3) = vals_phi_3;

            H_phi = sp_mat(Hphi_location, Hphi_values);
            HH_phi = H_phi.t() * H_phi;
        }

        // ------------------------------------------------------------
        // draw from condition for sigma^2_c
        // ------------------------------------------------------------
        S_hat_sigma2_c = S_sigma2_c + 0.5 * as_scalar((y - tau).t() * HH_phi * (y - tau));
        //S_hat_sigma2_c = 1 / S_hat_sigma2_c;
        sigma2_c = 1 / arma::randg<double>(distr_param(ny_hat_sigma2_c, 1/S_hat_sigma2_c));

        // ------------------------------------------------------------
        // draw from conditional for sigma^2_tau
        // ------------------------------------------------------------
        // gamma = (tau(0), tau(-1))'
        for (int j = 0; j < TT + 1; j++)
        {
            if (j > 1)
            {
                del_tau(j) = tau(j - 1) - tau(j - 2);
            }
            else if (j == 0)
            {
                del_tau(j) = gamma(0) - gamma(1);
            }
            else if (j == 1)
            {
                del_tau(j) = tau(j - 1) - gamma(0);
            }
        }
        sigma2_tau_grid = linspace(randu<double>() / 1000, b_sigma2_tau - randu<double>() / 1000, n_grid);
        lp_sigtau2 = -log(sigma2_tau_grid) * (TT / 2) - sum(square(diff(del_tau))) / (2 * sigma2_tau_grid);
        p_sigtau2 = exp(lp_sigtau2 - max(lp_sigtau2)) / sum(exp(lp_sigtau2 - max(lp_sigtau2)));
        cdf_sigtau2 = cumsum(p_sigtau2);
        idx = find(randu<double>() < cdf_sigtau2, 1);
        sigma2_tau = as_scalar(sigma2_tau_grid(idx));

        // ------------------------------------------------------------
        // draw from conditional for gamma
        // ------------------------------------------------------------
        K_gamma = invV_gamma + X_gamma.t() * (HH_2 * X_gamma) / sigma2_tau;
        XX_temp_gamma = invV_gamma * gamma_0 + X_gamma.t() * (HH_2 * tau) / sigma2_tau;
        gamma_hat = solve(K_gamma, XX_temp_gamma);
        gamma = gamma_hat + inv(chol(K_gamma, "upper")) * randn(2, 1);

        // ------------------------------------------------------------
        // Store draws after burn-in period
        // ------------------------------------------------------------
        if (i >= nburn)
        {
            tau_store.col(i - nburn) = tau;
            phi_store.col(i - nburn) = phi;
            gamma_store.col(i - nburn) = gamma;
            sigma2_c_store(i - nburn) = sigma2_c;
            sigma2_tau_store(i - nburn) = sigma2_tau;
        }
    }

    cout << "phi: " << mean(phi_store, 1) << endl;
    cout << "gamma: " << mean(gamma_store, 1) << endl;
    cout << "sigma2_tau: " << mean(sigma2_tau_store) << endl;
    cout << "sigma2_c: " << mean(sigma2_c_store) << endl;

    return 0;
}
