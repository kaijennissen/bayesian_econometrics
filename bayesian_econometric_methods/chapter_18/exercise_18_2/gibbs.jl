#!/usr/bin/env julia
#
# add packages
# 1) Pkg.add("PackageName")
# 2) type "]" to get into Package mode
# add PackageName
# leave Package mode by typing "ctrl + c"
using LinearAlgebra;
using CSV
using SparseArrays
#using Distributions
using Random
using Plots

function sparse_transpose(X::SparseMatrixCSC{Float64,Int64})
    i, j, v_ij = findnz(X);
    Xt = sparse(j, i, v_ij);
    return Xt
end

function sparse_diagonal(X::SparseMatrixCSC{Float64,Int64})
    i, j, v_ij = findnz(X);
    diag_elem = i.==j;
    new_i = i[diag_elem];
    new_v_ij = v_ij[diag_elem]
    Xt = sparse(new_i, new_i, new_v_ij);
    return Xt
end

function gibbs_sampler(y, nsim, burnin)

    T = length(y);

    # prior
    a0 = [750.0;750.0];
    B0 = 100*Diagonal([1,1]);
    phi0 = [1.3 -.7]';
    iVphi = Diagonal([1,1]);

    nu_sigc2 = 3.0;
    S_sigc2 = 2.;
    sigtau2_ub = .01;

    # initialize for storeage
    store_theta = zeros(nsim,6); # [phi, sigc2, sigtau2, tau0]
    store_tau = zeros(nsim,T);
    store_mu = zeros(nsim,T);    # annualized trend growth

    # initialize the Markov chain
    phi = [1.34 -.7]';
    tau0 = [y[1] y[1]]'; # [tau_{0}, tau_{-1}]
    sigc2 = .5;
    sigtau2 = .001;

    # construct a few things
    sparse_diag = sparse(1.0I, T, T)
    sparse_subdiag_1 = sparse(collect(2:T),collect(1:(T-1)),ones(T-1), T, T)
    sparse_subdiag_2 = sparse(collect(3:T),collect(1:(T-2)),ones(T-2), T, T)

    H2 = sparse_diag - 2 * sparse_subdiag_1 + sparse_subdiag_2;
    HH2 = H2' * H2;

    Hphi = sparse_diag - phi[1]*sparse_subdiag_1 - phi[2]*sparse_subdiag_2;
    HHphi = Hphi' * Hphi;

    Xtau0 = [collect(2:T+1) -collect(1:T)];
    n_grid = 500;
    count_phi = 0;

<<<<<<< Updated upstream
=======

>>>>>>> Stashed changes
    function f_tau(x, T, del_tau)
         return (-T/2*log.(x) - sum((del_tau[2:T] - del_tau[1:T-1]).^2)./(2*x));
    end

    for isim in 1:nsim+burnin
        # sample taup
        alp_tau = H2 \ [2*tau0[1]-tau0[2];-tau0[1];spzeros(T-2)];
        Ktau = triu(HH2 / sigtau2 + HHphi / sigc2);
        Ktau = Ktau + Ktau' - sparse_diagonal(Ktau);
        C = cholesky(Ktau, perm=1:size(Ktau, 1));
        tau_hat = Ktau \ (HH2 * alp_tau / sigtau2 + HHphi * y / sigc2);

        L = sparse(C.L);
        tau = tau_hat + L' \ rand(Normal(),T);

        # sample phi
        c = y - tau;
        Xphi = [[0;c[1:T-1]] [0;0;c[1:T-2]]];
        XXphi = Xphi' * Xphi;
        Kphi = iVphi + XXphi / sigc2;
        phi_hat = Kphi \ (iVphi * phi0 + Xphi' * c / sigc2);
        Cphi = cholesky(Kphi);
        phic = phi_hat + Cphi \ rand(Normal(),2);
        if sum(phic) < .99 && phic[2] - phic[1] < .99 && phic[2] > -.99
            phi = phic;
            Hphi = sparse_diag - phi[1]*sparse_subdiag_1 - phi[2]*sparse_subdiag_2;
            HHphi = Hphi'*Hphi;
            count_phi = count_phi + 1;
        end

        # sample sigc2
        G = Gamma(nu_sigc2 + T/2, 1/(S_sigc2 + .5*(c'*HHphi*c)[1]));
        sigc2 = 1 / rand(G);

        # sample sigtau2
        del_tau = [tau0[1];tau[1:T]] - [tau0[2];tau0[1];tau[1:T-1]];
        sigtau2_grid = collect(LinRange(rand(Uniform())/1000,sigtau2_ub-rand(Uniform())/1000,n_grid));
        lp_sigtau2 = f_tau(sigtau2_grid, T, del_tau);
        p_sigtau2 = exp.(lp_sigtau2 .- maximum(lp_sigtau2));
        p_sigtau2 = p_sigtau2/sum(p_sigtau2);
        cdf_sigtau2 = cumsum(p_sigtau2);
        sigtau2 = sigtau2_grid[rand(Uniform()).<cdf_sigtau2][1];

        # sample tau0
        Ktau0 = B0\Matrix(I,2 , 2) + Xtau0' * HH2 * Xtau0 / sigtau2;
        tau0_hat = Ktau0 \ (B0 \ a0 + Xtau0' * HH2 * tau / sigtau2);
        Ctau0 = cholesky(Ktau0);
        tau0 = tau0_hat + Ctau0.U \ randn(2);

        if isim > burnin
            i = isim - burnin;
            store_tau[i,:] = tau;
            store_theta[i,:] = [phi' sigc2 sigtau2 tau0'];
        end
    end

    tau_hat = mean(store_tau, dims = 1)';
    theta_hat = mean(store_theta, dims=1)';
    println(theta_hat);
    return store_tau, store_theta
end

data_raw = CSV.read("./usgdp.csv", header = 0, ignoreemptylines=true);
abc = reshape(Matrix(data_raw), 252);
data = 100*log.(abc);
y = data;

@time store_tau, store_theta  = gibbs_sampler(y, 500, 500);
@time store_tau, store_theta  = gibbs_sampler(y, 50000, 20000);
tau_hat = collect(mean(store_tau, dims = 1)');
theta_hat = mean(store_theta, dims=1)';

## plot of graphs
theta_CI = zeros(252,2)
for i in 1:252
    theta_CI[i,:] = collect(quantile(store_tau[:,i], (.025, .975)))
end

tt = (1947:.25:2009.75);
output_gap = reshape(y-tau_hat, 252);

output_gap = [y y y] - [theta_CI[:,1] tau_hat theta_CI[:,2]]
plot(tt, output_gap)
plot(output_gap[:,2], ribbon=(output_gap[:,2]-output_gap[:,1], output_gap[:,3]-output_gap[:,2]))
