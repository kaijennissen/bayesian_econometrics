#include <iostream>
#include <random>
#include <Eigen>
//#include <"Eigen/SparseCore">
//#include <"Eigen/SparseCholesky">


using namespace Eigen;
using namespace std;

int main()
{

    int nsim = 10;
    int nburn = 10;
    int TT = 252;
    int j;

    MatrixXd y = MatrixXd::Random(TT, 1) + VectorXd::LinSpaced(TT, 750, 750 + TT);

    // ----------------------------------------------------------------
    // Priors
    // ----------------------------------------------------------------
    // phi
    MatrixXd phi_0(2, 1);
    phi_0(0, 0) = 1.3;
    phi_0(1, 0) = -0.7;
    MatrixXd V_phi = MatrixXd::Identity(2, 2);
    MatrixXd invV_phi = V_phi.inverse();

    // gamma
    MatrixXd gamma_0 = MatrixXd::Constant(2, 1, 750);
    MatrixXd V_gamma = MatrixXd::Identity(2, 2) * 100;
    MatrixXd invV_gamma = V_gamma.inverse();

    // sigma2_tau
    double b_sigma2_tau = 0.01;

    // sigma2_c
    double ny_sigma2_c = 3;
    double S_sigma2_c = 2;

    // ----------------------------------------------------------------
    // Initial Values
    // ----------------------------------------------------------------
    MatrixXd phi(2, 1);
    phi(0, 0) = -1.34;
    phi(1, 0) = 0.7;
    double sigma2_c = {.5};
    double sigma2_tau = .001;
    MatrixXd tau_0 = MatrixXd::Constant(2, 1, y(0, 0));

    // H_2
    SparseMatrix<double> H_2(TT, TT);
    MatrixXd invH_2;
    for (j = 0; j < TT; j++)
    {
        H_2.insert(j, j) = 1;

        if (j > 0)
        {
            H_2.insert(j, j - 1) = -2;
        }
        if (j > 1)
        {
            H_2.insert(j, j - 2) = 1;
        }
    }
    H_2.makeCompressed();

    SparseMatrix<double> HH_2 = H_2.transpose() * H_2;

    SimplicialLDLT<SparseMatrix<double> > solverH2;
    solverH2.compute(H_2);
    //if(solverH2.info()!=Success) {
    // decomposition failed
    // return;
    //}
    // invH_2 = solverH2.solve(MatrixXd::Identity(TT));

    //if(solverH2.info()!=Success) {
    // solving failed
    //  return;
    //}

    // std::cout << "H_2=\n"
    //           << H_2 << std::endl;

    // H_phi
    SparseMatrix<double> H_phi(TT, TT);
    for (j = 0; j < TT; j++)
    {
        H_phi.insert(j, j) = 1;
        if (j > 0)
        {
            H_phi.insert(j, j - 1) = phi(0);
        }
        if (j > 1)
        {
            H_phi.insert(j, j - 2) = phi(1);
        }
    }
    H_phi.makeCompressed();
    SparseMatrix<double> HH_phi = H_phi.transpose() * H_phi;
    // std::cout << "H_phi=\n"
    //           << H_phi << std::endl;

    // X_gamma
    MatrixXd X_gamma(TT, 2);
    for (j = 0; j < TT; j++)
    {
        X_gamma(j, 0) = j + 2;
        X_gamma(j, 1) = -(j + 1);
    }
    // std::cout << "HX_gamma=\n"
    //           << X_gamma << std::endl;

    // ----------------------------------------------------------------
    // set up intermediate data structures and storage variables
    // ----------------------------------------------------------------

    // intermediate vars tau
    VectorXd S_sigma2_c_temp = MatrixXd::Zero(1, 1);
    MatrixXd alpha_tau_tilde = MatrixXd::Zero(TT, 1);
    MatrixXd alpha_tau(TT, 1);
    SparseMatrix<double> K_tau(TT, TT);
    MatrixXd XX_tau(TT, 1);
    MatrixXd Z_tau(TT, 1);
    MatrixXd tau_hat(TT, 1);
    MatrixXd tau(TT, 1);

    // intermediate vars phi
    MatrixXd K_phi(2, 2);
    MatrixXd XX(2, 2);
    MatrixXd phi_hat(2, 1);
    MatrixXd phic(2, 1);
    MatrixXd c(TT, 1);
    MatrixXd Z_phi(2, 1);
    MatrixXd X_phi(TT, 2);

    // sigma2_tau
    int n_grid = 500;
    MatrixXd del_tau(TT + 1, 1);
    MatrixXd sigma2_tau_grid(n_grid, 1);
    MatrixXd lp_sigtau2(n_grid, 1);
    MatrixXd p_sigtau2(n_grid, 1);
    MatrixXd cdf_sigtau2(n_grid, 1);
    //uvec idx;

    // gamma
    MatrixXd K_gamma(2, 2);
    MatrixXd XX_temp_gamma(2, 2);
    MatrixXd gamma_hat(2, 1);
    MatrixXd gamma(2, 1);

    // Store
    MatrixXd tau_store = MatrixXd::Zero(TT, nsim);
    MatrixXd phi_store = MatrixXd::Zero(TT, nsim);
    MatrixXd sigma2_c_store = MatrixXd::Zero(1, nsim);
    MatrixXd sigma2_tau_store = MatrixXd::Zero(1, nsim);
    MatrixXd gamma_store = MatrixXd::Zero(2, nsim);

    default_random_engine generator;
    normal_distribution<double> std_norm(0, 1);


    for (int i = 0; i < nsim+nburn; i++) {

    // ------------------------------------------------------------
    // draw from conditional for tau
    // ------------------------------------------------------------
    VectorXd v2;
    v2 = VectorXd::Zero(TT - 2);
    alpha_tau_tilde.col(0).head(2) << 2 * tau_0(1) - tau_0(0), -tau_0(1);
    alpha_tau = solverH2.solve(alpha_tau_tilde);
    K_tau = HH_phi / sigma2_c + HH_2 / sigma2_tau;

    SimplicialLLT<SparseMatrix<double> > solver_Ktau;
    solver_Ktau.compute(K_tau);
    MatrixXd Lt_K_tau = solver_Ktau.matrixU();

    XX_tau = HH_phi * y / sigma2_c + HH_2 * alpha_tau / sigma2_tau;
    tau_hat = solver_Ktau.solve(XX_tau);
    LLT<MatrixXd> llt_K_tau(K_tau);

    for (j = 0; j < TT; j++)
    {
        Z_tau(j, 0) = std_norm(generator);
    }
    tau = tau_hat + Lt_K_tau.inverse() * Z_tau;

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
    K_phi = invV_phi + X_phi.transpose() * X_phi / sigma2_c;
    XX = invV_phi * phi_0 + X_phi.transpose() * c / sigma2_c;

    phi_hat = K_phi.bdcSvd(ComputeThinU | ComputeThinV).solve(XX);
    MatrixXd Lt_K_phi = K_phi.llt().matrixU();

    // SimplicialLLT<SparseMatrix<double>> solver_K_phi;
    // solver_K_phi.compute(K_phi);
    // MatrixXd Lt_K_phi = solver_K_phi.matrixU();
    // phi_hat = solver_K_phi.solve(XX);

    for (j = 0; j < 2; j++)
    {
        Z_phi(j, 0) = std_norm(generator);
    }

    phic = phi_hat + Lt_K_phi.inverse() * Z_phi;

    if (phic(0, 0) + phic(1, 0) < .99 and phic(1, 0) - phic(0, 0) < .99 and phic(1, 0) > -.99)
    {
        phi = phic;
        // calculate H_phi in every iteration
        for (int j = 0; j < TT; j++)
        {
            H_phi.coeffRef(j, j) = 1;
            if (j > 0)
            {
                H_phi.coeffRef(j, j - 1) = phi(0);
            }
            if (j > 1)
            {
                H_phi.coeffRef(j, j - 2) = phi(1);
            }
        }
        HH_phi = H_phi.transpose() * H_phi;
    }

    // ------------------------------------------------------------
    // draw from condition for sigma^2_c
    // ------------------------------------------------------------
    S_sigma2_c_temp = 0.5 * (y - tau).transpose() * (H_phi * (y - tau));
    //ny_sigma2_c+TT/2 , 1/ (S_sigma2_c + S_sigma2_c_temp(0,0))
    gamma_distribution<double> gamma_dist(2,1);
    
    sigma2_c = 1/gamma_dist(generator);

    // std::cout << "sigma2_c=\n"
    //           << sigma2_c << std::endl;

    //         // ------------------------------------------------------------
    //         // draw from conditional for sigma^2_tau
    //         // ------------------------------------------------------------
    //         for (int j = 0; j<TT+1; j++){
    //             if (j>1){
    //                 del_tau(j,0) = tau(j-1)-tau(j-2);
    //             } else if (j==0){
    //                 del_tau(j,0) = tau_0(1)-tau_0(0);
    //             } else if(j==1){
    //                 del_tau(j,0) = tau(0)-tau_0(1);
    //             }
    //         }
    //         sigma2_tau_grid = linspace(as_scalar(randu(1)) / 100, b_sigma2_tau - as_scalar(randu(1)) / 100, n_grid);
    //         for (int j=0; j<n_grid; j++) {
    //             lp_sigtau2(j,0) = as_scalar(-log(sigma2_tau_grid(j)) * TT / 2 - sum(square(diff(del_tau))) / (2 * sigma2_tau_grid(j)));
    //         }
    //         p_sigtau2 = exp(lp_sigtau2 - as_scalar(max(lp_sigtau2))) / as_scalar(sum(exp(lp_sigtau2 - as_scalar(max(lp_sigtau2)))));
    //         cdf_sigtau2 = cumsum(p_sigtau2);
    //         idx = find(as_scalar(randu(1)) < cdf_sigtau2, 1);
    //         sigma2_tau = as_scalar(sigma2_tau_grid(idx));

    //         // ------------------------------------------------------------
    //         // draw from conditional for gamma
    //         // ------------------------------------------------------------
    //         mat K_gamma = invV_gamma + X_gamma.t() * (HH_2 * X_gamma) / as_scalar(sigma2_tau);
    //         mat XX_temp_gamma = invV_gamma * gamma_0 + X_gamma.t() * (HH_2 * tau) / as_scalar(sigma2_tau);
    //         mat gamma_hat = solve(K_gamma, XX_temp_gamma);
    //         mat gamma = gamma_hat + inv(chol(K_gamma)) * randn(2,1);

    //         // ------------------------------------------------------------
    //         // Store draws after burn-in period
    //         // ------------------------------------------------------------
    //         if (i >= nburn){
    //             tau_store.col(i-nburn) = tau;
    //             phi_store.col(i-nburn) = phi;
    //             sigma2_c_store.col(i-nburn) = sigma2_c;
    //             sigma2_tau_store.col(i-nburn) = sigma2_tau;
    //             gamma_store.col(i-nburn) = gamma;
    //             }

    //         }

    //     to_return(0) = tau_store;
    //     to_return(1) = phi_store;
    //     to_return(2) = sigma2_c_store;
    //     to_return(3) = sigma2_c_store;
    //     to_return(4) = gamma_store;

    //     return to_return;
    }
}
