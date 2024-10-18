"""
M-elektrodica:
                Kinetic functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import numpy as np
from .Constants import *

def k0 (Ea, DG, ne, beta, eta):
    return Ea - np.array([ne * beta * eta, DG + ne * (1 - beta) * eta])

def k (Ga, DG, ne, beta, V):
    return Ga - np.array([(1-beta)*(V*ne - DG),
                          -DG - beta * (V*ne - DG)])

def k3 (Ga, DG, ne, beta, V):
        return Ga - np.array([beta*(V*ne + DG),
                              DG + (1 - beta) * (V*ne + DG)])

def concentra (c_reactants, c_products, theta):
    return np.concatenate([c_reactants, c_products, theta, np.array([1-np.sum(theta[:-1])])])

def powerlaw (c: np.ndarray, upsilon: np.ndarray) -> np.ndarray:
    return np.array([np.prod(c ** (-upsilon * (upsilon < 0)), axis=1),
                    -np.prod(c ** (upsilon * (upsilon > 0)), axis=1)])

def rate (k: np.ndarray, r: np.ndarray) -> np.ndarray:
    return np.sum(k*r, axis=0)

def dc_dt (v: np.ndarray, upsilon: np.ndarray) -> np.ndarray:
    return np.dot(v, upsilon)

def aux(c_reactants: np.ndarray, c_products: np.ndarray, theta: np.ndarray, V, Ga, DG, beta, ne, nu, T):
    kexp = np.exp(-k(Ga, DG, ne, beta, V) / k_B / T)
    c = concentra(c_reactants, c_products, theta)
    r = powerlaw(c, nu)
    return rate(kexp, r)

def current(c_reactants: np.ndarray, c_products: np.ndarray, theta: np.ndarray, V, Ga, DG, beta, ne, nu, T):
    return np.dot(ne, aux(c_reactants, c_products, theta, V, Ga, DG, beta, ne, nu, T)) * F * k_B * T / h






