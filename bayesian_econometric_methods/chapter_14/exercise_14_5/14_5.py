# This m-file estimates an ordered probit
# model using generated data, and the
# algorithm of Albert and Chib(1993).
import numpy as np
import random
from scipy.stats import truncnorm, norm
from copy import deepcopy
import matplotlib.pyplot as plt
random.seed(42)

# Set the parameters of the generated data
nobs: int = 1500
L: int = 3
beta_0: float = .25
beta_1: float = 1.5
x: np.ndarray = np.random.randn(nobs, 1)
z: np.ndarray = beta_0 + beta_1 * x + np.random.randn(nobs, 1)
big_x: np.ndarray = np.append(np.ones((nobs, 1)), x, axis=1)

alpha_2: float = .8
alpha_3: float = 1.1
# create the observed y variable based on the
# latent data
y: np.ndarray = np.zeros((nobs, 1))
D1: np.ndarray = np.zeros((nobs, 1))
D2: np.ndarray = np.zeros((nobs, 1))
D3: np.ndarray = np.zeros((nobs, 1))
D4: np.ndarray = np.zeros((nobs, 1))

for i, zz in enumerate(z):
    if (zz < 0):
        y[i, 0] = 1
        D1[i, 0] = 1
    elif (zz > 0) & (zz < alpha_2):
        y[i, 0] = 2
        D2[i, 0] = 1
    elif (zz > alpha_2) & (zz < alpha_3):
        y[i, 0] = 3
        D3[i, 0] = 1
    else:
        y[i, 0] = 4
        D4[i, 0] = 1

points_2 = np.nonzero(D2)[0]
points_3 = np.nonzero(D3)[0]
points_4 = np.nonzero(D4)[0]

# intial conditions
alph1: float = 1
alph2: float = 1.5
betas: np.ndarray = np.zeros([2, 1])
n_iter: int = 7500
burn: int = 2500
betas_final = np.zeros((n_iter - burn, big_x.shape[1]))
alph_final = np.zeros((n_iter - burn, 2))


def trunc_normal(mu, var, direct):

    n = mu.shape[0]
    uniforms = np.random.uniform(size=(n, 1))
    stderrs = np.sqrt(var)

    c = norm.cdf(-mu / stderrs)
    p = (c * (1 - direct) + (1 - c) * direct) * uniforms + c * direct
    draws = mu + stderrs * norm.ppf(p)

    return draws


def trunc2_normal(mu, var, a, b):

    n = mu.shape[0]
    stderrs = np.sqrt(var)

    if a == -999:
        a_term = np.zeros((n, 1))
    else:
        a_term = norm.cdf((a - mu) / stderrs)

    if b == 999:
        b_term = np.ones((n, 1))
    else:
        b_term = norm.cdf((b - mu) / stderrs)

    uniforms = np.random.uniform(size=(n, 1))

    p = a_term + uniforms * (b_term - a_term)
    draws = mu + stderrs * norm.ppf(p)

    return draws


for i in range(n_iter):

    # posterior conditional for z
    zaug = (D1 * trunc_normal(mu=big_x @ betas, var=1, direct=0)
            + D2 * trunc2_normal(mu=big_x @ betas, var=1, a=0, b=alph1)
            + D3 * trunc2_normal(mu=big_x @ betas, var=1, a=alph1, b=alph2)
            + D4 * trunc2_normal(mu=big_x @ betas, var=1, a=alph2, b=999))

    # posterior conditional for beta
    D_beta = np.linalg.inv(big_x.T @ big_x)
    d_beta = big_x.T @ zaug
    H_beta = np.linalg.cholesky(D_beta)
    betas = D_beta @ d_beta + H_beta.T @ np.random.randn(big_x.shape[1], 1)

    # posterior conditional for alpha
    # alpha1
    m1 = np.max([0, np.max(zaug[points_2])])
    m2 = np.min(zaug[points_3])
    u = np.random.uniform(low=0., high=1.)
    alph = m1 + (m2 - m1) * u

    # alpha2
    m1 = np.max([0, np.max(zaug[points_3])])
    m2 = np.min(zaug[points_4])
    u = np.random.uniform(low=0., high=1.)
    alph2 = m1 + (m2 - m1) * u

    if (i > burn):
        betas_final[i - burn, :] = betas.T
        alph_final[i - burn, :] = np.array([alph1, alph2])


print(f"betas: {np.mean(betas_final, axis=0)}")
print(f"alphas: {np.mean(alph_final, axis=0)}")
