import numpy as np
from scipy.stats import poisson

def bell_pmf(k, mu):
    return poisson.pmf(k, mu)

def bell_rvs(mu, size=1):
    return poisson.rvs(mu, size=size)
