import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from .bell import bell_pmf

def g_density(y, mu, w):
    return np.where(y == 0, w + (1 - w) * bell_pmf(y, mu), (1 - w) * bell_pmf(y, mu))

def ell_alpha(y_d, mu, w, alpha, N=20):
    y_vals = np.arange(0, N + 1)
    S = np.sum(g_density(y_vals, mu, w) ** (1 + alpha))
    return S - (1 + 1 / alpha) * (g_density(y_d, mu, w) ** alpha)

def fonc_H(theta, Y, X, Z, alpha):
    b = theta[:2]
    g = theta[2]
    w_i = expit(g * Z)
    mu_i = np.exp(b @ X)
    values = np.array([ell_alpha(Y[i], mu_i[i], w_i[i], alpha) for i in range(len(Y))])
    return np.mean(values)

def estime_rep(Y, X, Z, alpha):
    Theta_In = np.array([0.1, 0.1, 0.1])
    lik_opt = lambda theta: fonc_H(theta, Y, X, Z, alpha)
    result = minimize(lik_opt, Theta_In, method='L-BFGS-B')
    return result.x

def monte_carlo(donnees, alpha, true_b, true_g):
    True_theta = np.concatenate([true_b, [true_g]])
    My = donnees['My']
    Mx = donnees['Mx']
    Z = np.ones(My.shape[1])
    Mat_Est = np.zeros((My.shape[0], len(True_theta)))
    for k in range(My.shape[0]):
        print(f'Estimation réplication {k+1}/{My.shape[0]}')
        Mat_Est[k, :] = estime_rep(My[k, :], Mx[:, :, k], Z, alpha)
    Mean = np.mean(Mat_Est, axis=0)
    biais = Mean - True_theta
    rmse = np.sqrt(np.mean((Mat_Est - True_theta) ** 2, axis=0))
    return {'True_theta': True_theta, 'Mean.est': np.round(Mean, 3), 'biais': np.round(biais, 3), 'Rmse': np.round(rmse, 3)}
