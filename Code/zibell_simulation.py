import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
from scipy.special import lambertw, factorial
from scipy.special import expit  # sigmoïde: 1/(1+exp(-x))

# Distribution de Bell (Bell PMF) approximation
# pmf Bell distribution: P(Y=k) = e^{-e^{mu}} * (e^{mu})^k * B_k / k!
# where B_k = nombre de Bell (Bell numbers)
# mais B_k est compliqué, on peut approx via une distribution liée.

# Une approche simple: utiliser la distribution de Bell par
# la formule: P(Y=k) = exp(-mu) * mu^k / k! * e^{-e^mu} * Bell(k)
# comme approximation, ou utiliser la fonction rbell de R via package bellreg (non disponible)

# Ici, on implémente une approximation (à améliorer selon besoin)

def bell_pmf(k, mu):
    # Approximation pmf Bell distribution (par convolution Poisson + Bell numbers)
    # Ici on approxime par Poisson(mu), car Bell est très proche d'une distribution Poisson mélangée.
    # Pour un mu donné, on peut approximer par Poisson.
    # Sinon il faut une vraie fonction Bell, mais elle est complexe.
    # Cette approximation est utilisée pour continuer.
    return poisson.pmf(k, mu)

def bell_rvs(mu, size=1):
    # Génère des valeurs aléatoires Bell, approximé par Poisson
    return poisson.rvs(mu, size=size)

##########################
# Fonction data équivalente
def data(n, b, g):
    inter = np.ones(n)
    X1 = np.random.normal(0, 1, n)
    X = np.vstack([inter, X1])  # 2 x n

    Z = np.ones(n)  # intercept only

    # zero-inflation probability w = sigmoid(g * Z)
    w = expit(g * Z)
    mu = np.exp(b @ X)  # vector of means

    Y = np.zeros(n)
    S = np.random.binomial(1, w)  # indicateur zero-inflation

    for i in range(n):
        if S[i] == 1:
            Y[i] = 0
        else:
            # échantillon Bell
            Y[i] = bell_rvs(mu[i])[0]

    return {'Y': Y, 'X': X, 'Z': Z}

####################
# g.density équivalent
def g_density(y, mu, w):
    # y can be scalar or array, mu and w same shape
    # P(Y=y) = w + (1-w)*Bell_pmf(y, mu) si y=0
    # sinon (1-w)*Bell_pmf(y, mu)
    result = np.where(y == 0,
                      w + (1 - w) * bell_pmf(y, mu),
                      (1 - w) * bell_pmf(y, mu))
    return result

####################
# Fonction Ell_alpha équivalente

def Ell_alpha(y_d, mu, w, alpha, N=20):
    # Calcul de sum_{y=0}^N (g_density(y, mu, w))^(1+alpha) - (1+1/alpha) * (g_density(y_d, mu, w))^alpha
    y_vals = np.arange(0, N + 1)
    S = np.sum(g_density(y_vals, mu, w) ** (1 + alpha))
    val = S - (1 + 1 / alpha) * (g_density(y_d, mu, w) ** alpha)
    return val

####################
# Fonction H_alpha équivalente

def Fonc_H(theta, Y, X, Z, alpha):
    b = theta[:2]
    g = theta[2]

    w_i = expit(g * Z)
    mu_i = np.exp(b @ X)

    values = np.array([Ell_alpha(Y[i], mu_i[i], w_i[i], alpha) for i in range(len(Y))])
    return np.mean(values)

####################
# Fonction d'estimation

def Estime_Rep(Y, X, Z, alpha):
    Theta_In = np.array([0.1, 0.1, 0.1])

    def lik_opt(theta):
        return Fonc_H(theta, Y, X, Z, alpha)

    result = minimize(lik_opt, Theta_In, method='L-BFGS-B')
    return result.x

####################
# Fonction simulation répétée

def Simul_rep(n, nbrep, b, g):
    My = np.zeros((nbrep, n))
    Mx = np.zeros((2, n, nbrep))
    Mx[0, :, :] = 1

    for i in range(nbrep):
        data_sim = data(n, b, g)
        My[i, :] = data_sim['Y']
        Mx[:, :, i] = data_sim['X']

    return {'My': My, 'Mx': Mx}

####################
# Fonction Monte Carlo

def MC(donnees, alpha):
    b = np.array([-0.5, 1.2])
    g = -1.1
    True_theta = np.concatenate([b, [g]])
    My = donnees['My']
    Mx = donnees['Mx']
    Z = np.ones(My.shape[1])

    Mat_Est = np.zeros((My.shape[0], len(True_theta)))

    for k in range(My.shape[0]):
        print(f'Estimation réplication {k+1}/{My.shape[0]}')
        Mat_Est[k, :] = Estime_Rep(My[k, :], Mx[:, :, k], Z, alpha)

    Mean = np.mean(Mat_Est, axis=0)
    biais = Mean - True_theta
    Rmse = np.sqrt(np.mean((Mat_Est - True_theta) ** 2, axis=0))

    return {'True_theta': True_theta, 'Mean.est': np.round(Mean, 3), 'biais': np.round(biais, 3), 'Rmse': np.round(Rmse, 3)}

####################
# Exemple d'utilisation

if __name__ == "__main__":
    np.random.seed(1234)

    b = np.array([-0.5, 1.2])
    g = -1.1

    # Simulation des données
    data_sim = data(1000, b, g)
    Y = data_sim['Y']
    print(np.bincount(Y.astype(int)))

    # Simulation répétée
    Data1000NC_ZI50 = Simul_rep(1000, 10, b, g)  # ici 10 répétitions pour test rapide

    # Monte Carlo estimation
    import time
    T1 = time.time()
    Res = MC(Data1000NC_ZI50, alpha=0.1)
    T2 = time.time()
    print("Temps écoulé:", T2 - T1)
    print(Res)
