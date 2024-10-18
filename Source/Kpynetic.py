"""
M-elektrodica:
                Kinetic functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import numpy as np
from .Constants import *

def k (Ga, DG, ne, beta, V):
    return - Ga + np.array([-beta, 1-beta]) * (V*ne+DG)

def concentra (c_reactants, c_products, theta):
    return np.concatenate([c_reactants, c_products, theta, np.array([1-np.sum(theta[:-1])])])

def powerlaw (c: np.ndarray, upsilon: np.ndarray) -> np.ndarray:
    return np.array([np.prod(c ** (-upsilon * (upsilon < 0)), axis=1),
                    -np.prod(c ** (upsilon * (upsilon > 0)), axis=1)])

def rate (k: np.ndarray, r: np.ndarray) -> np.ndarray:
    return np.sum(k*r, axis=0)

def dc_dt (v: np.ndarray, upsilon: np.ndarray) -> np.ndarray:
    return np.dot(v, upsilon)










