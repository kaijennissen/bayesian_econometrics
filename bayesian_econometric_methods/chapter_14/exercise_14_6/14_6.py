# This m-file estimates an ordered probit
# model using generated data, and the
# algorithm of Albert and Chib(1993).
import numpy as np
import random
from scipy.stats import truncnorm, norm, invgamma
from copy import deepcopy
import matplotlib.pyplot as plt
random.seed(42)


# Set the parameters of the generated data
nobs: int = 750
L: int = 3
beta_0: float = .3
beta_1: float = .85
x: np.ndarray = np.random.randn(nobs, 3)
z: np.ndarray = beta_0 + beta_1 * x + np.random.randn(nobs, 1)
big_x: np.ndarray = np.append(np.ones((nobs, 1)), x, axis=1)

alpha_2: float = 0.65
# create the observed y variable based on the
y: np.ndarray = np.zeros((nobs, 1))
D1: np.ndarray = np.zeros((nobs, 1))
D2: np.ndarray = np.zeros((nobs, 1))
D3: np.ndarray = np.zeros((nobs, 1))

for i, zz in enumerate(z):
    if (zz < 0):
        y[i, 0] = 1
        D1[i, 0] = 1
    elif (zz > 0) & (zz < alpha_2):
        y[i, 0] = 2
        D2[i, 0] = 1
    else:
        y[i, 0] = 3.
        D3[i, 0] = 1

points_2 = np.nonzero(D2)[0]
points_3 = np.nonzero(D3)[0]


# intial conditions
delt: float = 1.
alph: float = .5
betas: np.ndarray = np.zeros([2, 1])
n_iter: int = 5000
burn: int = 1500
betas_final = np.zeros((n_iter - burn, big_x.shape[1]))
alph_final = np.zeros((n_iter - burn, 1))

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
    zaug = (D1 * trunc_normal(mu=big_x @ betas, var=delt, direct=0)
            + D2 * trunc2_normal(mu=big_x @ betas, var=delt, a=0., b=1.)
            + D3 * trunc2_normal(mu=big_x @ betas, var=delt, a=1., b=999.))

    # posterior conditional for beta
    D_beta = np.linalg.inv((big_x.T @ big_x) / delt)
    d_beta = big_x.T @ zaug / delt
    H_beta = np.linalg.cholesky(D_beta)
    betas = D_beta @ d_beta + H_beta.T @ np.random.randn(big_x.shape[1], 1)

    # posterior conditional for delta
    shape = (nobs + big_x.shape[1]) / 2
    scale = (.5 * (zaug - big_x @ betas).T @ (zaug - big_x @ betas))[0, 0]
    delt = invgamma.rvs(a=shape, scale=scale)

    alph = 1 / np.sqrt(delt)
    betas_use = betas * alph

    if (i > burn):
        betas_final[i - burn, :] = betas_use.T
        alph_final[i - burn, 0] = alph


# posterior predictions


print(f"beta: {np.mean(betas_final, axis=0)}")
print(f"alpha: {np.mean(alph_final, axis=0)}")
