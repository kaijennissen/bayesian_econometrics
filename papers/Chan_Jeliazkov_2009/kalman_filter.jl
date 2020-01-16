#!/usr/bin/env julia
using LinearAlgebra
using CSV
using SparseArrays
using Distributions
using Random
using Plots
#using BandedMatrices
#using Profile


function sparse_transpose(X)
    i, j, v_ij = findnz(X)
    Xt = sparse(j, i, v_ij)
    return Xt
end

function SUR(X)
    r, c = size(X);
    idi = reshape((kron((1:r), ones(c, 1) )), r*c);
    idj = collect(1:r*c);
    X_SUR = sparse(idi, idj, reshape(X',r*c));
    return X_SUR
end

nsim = 10;
burnin = 10;

data_raw = CSV.read("USdata.csv", header = 0);
y = Matrix(data_raw);

function gibbs_sampler(y, nsim, burnin)

    y0 = y[1:2, :];
    Y = y[4:end, :];

    tt = size(y, 1);
    nn = size(y, 2);
    TT = size(Y, 1);
    qq = nn*(nn+1);
    TTqq = TT*qq;

    Y = reshape(Y',:,1);

    # priors #-----------------------------------------------------------------
    # Omega_11
    nu01 = nn + 3. # nn?
    S01 = I(nn)

    # Omega_22
    DD = 5 * I(qq);
    DD_inv = inv(DD);
    nu02 = 6 * ones(qq);
    S02 = 0.01 * ones(qq);

    # initial values #---------------------------------------------------------
    new_nu01 = nu01 + TT;
    new_nu02 = broadcast(+ , nu02/2, (TT-1)/2)[:];

    # Omega11
    Omega11 = cov(y);
    Omega11_inv = inv(Omega11);

    # H
    H1 = I(TTqq);
    H2 = [[zeros(qq, TTqq-qq); I(TTqq-qq)] zeros(TTqq, qq)];
    H = H1 - H2;
    HT = sparse_transpose(H)

    # S
    Omega22 = 0.01*I(qq);
    Omega22_inv = inv(Omega22);
    S = kron(I(TT), Omega22);
    S[1:qq,1:qq] = DD;

    S_inv = kron(I(TT), Omega22_inv);
    S_inv[1:qq,1:qq] = DD_inv;

    # G
    G = SUR([ones(TT*nn,1) kron(y[3:end-1, :], ones(nn))]);
    GT = sparse_transpose(G);

    # initialize for storeage
    store_beta = zeros(nsim, TTqq);
    store_Omega11 = zeros(nn, nn,  nsim);
    store_Omega22 = zeros(nsim, qq);

    for isim in 1:nsim+burnin

        S_inv = kron(sparse(I,TT,TT), Omega22_inv);
        S_inv[1:qq,1:qq] = DD_inv;
        K = HT * S_inv * H;

        #GT_Omega11_inv = GT * kron(sparse(I, TT, TT), Omega11_inv);
        GGL = tril(GT * kron(sparse(I, TT, TT), Omega11_inv) * G);
        GT_Omega11_inv_G = GGL + sparse_transpose(GGL) - Diagonal(GGL);
        GT_Omega11_inv_Y = GT * (kron(sparse(I, TT, TT), Omega11_inv) * Y);
        P = K + GT_Omega11_inv_G;

        # cholesky 1
        C = cholesky(P, perm=1:TTqq);
        L = sparse(C.L);
        #L * L' ≈ P;

        # Ctau = cholesky(P);
        # Ltau = sparse(Ctau.L);
        # #nrows_tau = size(Ctau, 1);
        # Ptau = sparse(1:TTqq, Ctau.p, ones(TTqq));
        # LtPtau = Ltau' * Ptau;
        # LtPtauT = sparse_transpose(LtPtau);

        # cholesky 3
        # C = cholesky(Matrix(P));
        # L = C.L;
        # L*L' ≈  Matrix(P);


        # beta
        # fwd = fwdTriSolve(LtPtau, GT_Omega11_inv_Y);
        # bwd = bwdTriSolve(LtPtauT, fwd);

        beta_hat = L'\(L\GT_Omega11_inv_Y);
        beta = beta_hat + L' \ rand(Normal(),TTqq);

        # Omega11
        e1 = reshape(Y - G * beta, nn,:);
        new_S01 = S01 + e1*e1';
        Omega11 = rand(InverseWishart(new_nu01, new_S01));
        Omega11_inv = sparse(Symmetric(Omega11\I(nn)));

        # Omega22
        e2 = reshape(H * beta, qq, TT)';
        new_S02 = (S02 + sum(e2[2:end,:].^2, dims=1)[:])/2;
        for i in 1:qq
            Omega22[i,i] = rand(InverseGamma(new_nu02[i], new_S02[i]));
        end
        Omega22_inv = Omega22\I(qq);

        # store
        if isim > burnin
            i = isim - burnin;
            store_beta[i,:] = beta;
            store_Omega11[:,:, i] = Omega11;
            store_Omega22[i,:] = diag(Omega22);
        end
    end

    return store_beta, store_Omega11, store_Omega22
end

# run first with a small number for compiltaion
@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 100, 100);
@time store_beta, store_Omega11, store_Omega22  = gibbs_sampler(y, 50000, 5000);2.54/200

beta_hat = mean(store_beta, dims = 1)';
Omega11_hat = mean(store_Omega11, dims = 3)[:,:,1];
Omega22_hat = mean(store_Omega22, dims = 1)';

beta = reshape(beta_hat, 4, 5, :);

l = @layout [a b ; c d]
p1 = plot(1:245, beta[1,:,:]', legend=false)
p2 = plot(1:245, beta[2,:,:]', legend=false)
p3 = plot(1:245, beta[3,:,:]', legend=false)
p4 = plot(1:245, beta[4,:,:]', legend=false)
plot(p1, p2, p3, p4, layout = l)


savefig( "fg")
