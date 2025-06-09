import numpy as np
from scipy.special import expit
from .bell import bell_rvs

def data(n, b, g):
    inter = np.ones(n)
    X1 = np.random.normal(0, 1, n)
    X = np.vstack([inter, X1])
    Z = np.ones(n)
    w = expit(g * Z)
    mu = np.exp(b @ X)
    Y = np.zeros(n)
    S = np.random.binomial(1, w)
    for i in range(n):
        Y[i] = 0 if S[i] == 1 else bell_rvs(mu[i])[0]
    return {'Y': Y, 'X': X, 'Z': Z}

def simul_rep(n, nbrep, b, g):
    My = np.zeros((nbrep, n))
    Mx = np.zeros((2, n, nbrep))
    Mx[0, :, :] = 1
    for i in range(nbrep):
        data_sim = data(n, b, g)
        My[i, :] = data_sim['Y']
        Mx[:, :, i] = data_sim['X']
    return {'My': My, 'Mx': Mx}
