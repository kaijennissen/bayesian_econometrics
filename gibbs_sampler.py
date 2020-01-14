import tensorflow as tf
import tensorflow_probability as tfp
import pandas as pd
import numpy as np
import scipy as scp
import scipy.stats as scps
import scipy.linalg as slin
import scipy.sparse as ssp
import scipy.sparse.linalg as ssplin
import timeit

tfd = tfp.distributions

nsim = 50
nburn = 10
total_runs = nsim + nburn
data = pd.read_csv("./papers/Chan_Jeliazkov_2009/USdata.csv", header=None)

y = data.values




# initial values #-------------------------------------------------------------

# Omega_11
Omega11 = np.cov(y.T)  # tfp.stats.covariance(y)
Omega11_inv = np.linalg.inv(Omega11)  # tf.linalg.inv(Omega11)
# H1 = tf.sparse.eye(
#     TT*qq,
#     num_columns=None,
#     dtype=tf.dtypes.float32,
# )
# H2 = tf.concat([
#     tf.zeros((qq,TT*qq)),
#     tf.sparse.eye(
#         num_rows=(TT-1)*qq,
#         num_columns=TT*qq,
#         dtype=tf.dtypes.float32)], axis=0)
# H = H1 - H2
# TODO tf.linalg.diag is fixed in nightly build
H1 = ssp.spdiags(np.ones(TT*qq), diags=0, m=TT*qq, n=TT*qq)  # tf.linalg.diag(np.ones(TT*qq), k=0, num_rows=TT+qq, num_cols=TT*qq)
H2 = ssp.spdiags(np.ones((TT-1)*qq), diags=-qq, m=TT*qq, n=TT*qq)  # tf.linalg.diag(np.ones((TT-1)*qq), k=-qq, num_rows=TT*qq, num_cols=TT*qq)
H = H1+H2
# H = tf.sparse.from_dense(H.toarray())  # tf.sparse(H)

# S
DD_inv = np.linalg.solve(DD, np.eye(DD.shape[0]))
Omega22 = 0.01 * np.eye(qq) # initial values for variance of betas?
Omega22_inv = np.linalg.solve(Omega22, np.eye(Omega22.shape[0]))
Omega22_inv = ssp.coo_matrix(Omega22_inv)

args_S = [DD] + [Omega22 for i in range(TT-1)]
S = slin.block_diag(*args_S)
S = ssp.coo_matrix(S)
# S = tf.sparse.from_dense(S.toarray())

args_S_inv = [DD_inv] + [Omega22_inv.toarray() for i in range(TT-1)]
S_inv = slin.block_diag(*args_S_inv)
S_inv = ssp.coo_matrix(S_inv)
# S_inv = tf.sparse.from_dense(S_inv.toarray())

# G
args_G = (np.kron(np.eye(nn), np.append([1], y[i,:])) for i in range(3, tt))
G = slin.block_diag(*args_G)
G = ssp.coo_matrix(G)
# G = tf.sparse.from_dense(G)

# store
store_beta = np.empty((nsim, TT*qq))
store_Omega11 = np.empty((nn, nn, nsim))
store_Omega22 = np.empty((nsim, qq))

class Gibbs_Sampler():

    def __init__(self, y, nsim, nburn):
        self._y = y
        self._nsim = nsim
        self._nbunr = nburn

        # set up
        y0 = y[0:3, ]
        Y = y[3:, :].ravel().reshape((-1, 1))
        self._tt = y.shape[0]
        self._nn = y.shape[1]
        self._TT = int(Y.shape[0] / self._nn)
        self._qq = self._nn * (self._nn + 1)

        # priors #---------------------------------------------------------------------
        # Omega_11
        self.nu01 = self._nn + 3.  # tf.constant(nn + 3.)  #
        self.S01 = np.eye(self._nn)  # tf.eye(2)

        # Omega_22
        self.DD = 5. * np.eye(self._qq)  # tf.eye(qq, dtype=tf.dtypes.float32)
        self.nu02 = 6. * np.ones(self._qq)  # tf.ones(qq, dtype=tf.dtypes.float32)
        self.S02 = 0.01 * np.ones(self._qq)  # tf.ones(qq, dtype=tf.dtypes.float32)

        # initial values #-------------------------------------------------------------
        # Omega_11
        Omega11 = np.cov(y.T)  # tfp.stats.covariance(y)
        Omega11_inv = np.linalg.inv(Omega11)  # tf.linalg.inv(Omega11)
        # H1 = tf.sparse.eye(
        #     TT*qq,
        #     num_columns=None,
        #     dtype=tf.dtypes.float32,
        # )
        # H2 = tf.concat([
        #     tf.zeros((qq,TT*qq)),
        #     tf.sparse.eye(
        #         num_rows=(TT-1)*qq,
        #         num_columns=TT*qq,
        #         dtype=tf.dtypes.float32)], axis=0)
        # H = H1 - H2
        # TODO tf.linalg.diag is fixed in nightly build
        H1 = ssp.spdiags(np.ones(TT * qq), diags=0, m=TT * qq,
                         n=TT * qq)  # tf.linalg.diag(np.ones(TT*qq), k=0, num_rows=TT+qq, num_cols=TT*qq)
        H2 = ssp.spdiags(np.ones((TT - 1) * qq), diags=-qq, m=TT * qq,
                         n=TT * qq)  # tf.linalg.diag(np.ones((TT-1)*qq), k=-qq, num_rows=TT*qq, num_cols=TT*qq)
        H = H1 + H2
        # H = tf.sparse.from_dense(H.toarray())  # tf.sparse(H)

        # S
        DD_inv = np.linalg.solve(DD, np.eye(DD.shape[0]))
        Omega22 = 0.01 * np.eye(qq)  # initial values for variance of betas?
        Omega22_inv = np.linalg.solve(Omega22, np.eye(Omega22.shape[0]))
        Omega22_inv = ssp.coo_matrix(Omega22_inv)

        args_S = [DD] + [Omega22 for i in range(TT - 1)]
        S = slin.block_diag(*args_S)
        S = ssp.coo_matrix(S)
        # S = tf.sparse.from_dense(S.toarray())

        args_S_inv = [DD_inv] + [Omega22_inv.toarray() for i in range(TT - 1)]
        S_inv = slin.block_diag(*args_S_inv)
        S_inv = ssp.coo_matrix(S_inv)
        # S_inv = tf.sparse.from_dense(S_inv.toarray())

        # G
        args_G = (np.kron(np.eye(nn), np.append([1], y[i, :])) for i in range(3, tt))
        G = slin.block_diag(*args_G)
        G = ssp.coo_matrix(G)
        # G = tf.sparse.from_dense(G)

        # store
        store_beta = np.empty((nsim, TT * qq))
        store_Omega11 = np.empty((nn, nn, nsim))
        store_Omega22 = np.empty((nsim, qq))


    def draw_beta(self:
        # beta --------------------------------------------------------------------
        args_S_inv = [DD_inv] + [Omega22_inv.toarray() for i in range(TT-1)]
        S_inv = slin.block_diag(*args_S_inv)
        S_inv = ssp.coo_matrix(S_inv)
        #S_inv = tf.sparse.from_dense(S_inv.toarray())
        #S_inv = tf.sparse.to_dense(S_inv)

        #H = tf.sparse.to_dense(H)
        #HT = tf.transpose(H)
        K = H.T * S_inv * H  # tf.linalg.matmul(HT, tf.linalg.matmul(S_inv, H))
        G_Omega11_inv = G.T * ssp.kron(ssp.eye(TT), Omega11_inv)
        G_Omega11_inv_G = G_Omega11_inv * G
        P = K + G_Omega11_inv_G
        L = ssp.csr_matrix(np.linalg.cholesky(P.toarray()))

        beta_hat = ssplin.spsolve_triangular(L.T, ssplin.spsolve_triangular(L,(G_Omega11_inv*Y), lower=True), lower=False)
        beta = beta_hat + ssplin.spsolve_triangular(L, np.random.normal(size=(beta_hat.shape[0],1)), lower=True)

        return beta


    def draw_Omega11():
        e1 = np.array(Y - G * beta).reshape((-1, nn)).T
        new_nu1 = nu1 + TT
        new_S1 = S1 + np.dot(e1, e1.T)

        Omega11 = scps.invwishart(df=new_nu1, scale=new_S1).rvs(size=1)
        Omega11_inv = np.linalg.inv(Omega11)

        return Omega11


    def draw_Omega22():
        e2 = np.array(H.dot(beta)).reshape((-1, TT)).T
        new_nu2 = (nu2 + (TT-1)) / 2
        new_S2 = (S2 + np.sum(np.square(e2)[1:, :], axis=0, keepdims=True))/2
        np.fill_diagonal(Omega22, scps.invgamma(a=new_nu2, scale=new_S2).rvs(size=20))

        return Omega22



timeit.timeit(draw_beta, number=20)/20
timeit.timeit(draw_Omega11, number=20)/20
timeit.timeit(draw_Omega22, number=20)/20