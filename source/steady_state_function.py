"""
M-elektrodica:
                Steady state functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""
import numpy as np
from scipy.optimize import fsolve
from kinetic_functions import dc_dt, aux, current

def object(variables: np.ndarray, V: float,
           Ga: np.ndarray, DG: np.ndarray, beta: np.ndarray,
           ne: np.ndarray, nu: np.ndarray, T: float,
           c0_reactants: np.ndarray, c0_products: np.ndarray,
           Fv: float, A: float):

     c_reactants = variables[:len(c0_reactants)]
     c_products  = variables[len(c0_reactants) : len(c0_reactants) + len(c0_products)]
     theta       = variables[len(c0_reactants) + len(c0_products):]


     rhs = np.concatenate([(c_reactants - c0_reactants)*Fv/A,
                           (c_products  - c0_products)*Fv/A,
                           np.zeros(len(theta))])

     return dc_dt(aux(c_reactants, c_products, theta, V, Ga, DG, beta, ne, nu, T), nu) - rhs


def steady_state(c0_reactants: np.ndarray, c0_products: np.ndarray, adsorbed: np.ndarray,
                 E: np.ndarray, Ga: float, DG: float, beta: float, ne: float,
                 nu: np.ndarray, T: float, FV: float, A: float) -> tuple:

    lreactants: int = len(c0_reactants)
    lproducts: int = len(c0_products)
    ladsorbed: int = len(adsorbed)
    lE: int = len(E)

    c_reactants0: np.ndarray = np.ones(lreactants)
    c_products0: np.ndarray = np.zeros(lproducts)
    theta0: np.ndarray = np.zeros(ladsorbed)

    init: np.ndarray = np.concatenate([c_reactants0, c_products0, theta0])


    c_reactants: np.ndarray = np.zeros((lE, lreactants))
    c_products: np.ndarray = np.zeros((lE, lproducts))
    theta: np.ndarray = np.zeros((lE, ladsorbed))
    fval: np.ndarray = np.zeros((lE, lreactants + lproducts + ladsorbed))
    j: np.ndarray = np.zeros(lE)

    for i in range(lE):
        steady_state_func = lambda  variables: object( variables, E[i],
                                                      Ga, DG, beta,
                                                      ne, nu[:(lreactants + lproducts + ladsorbed)], T,
                                                      c0_reactants, c0_products, FV, A)
        sol = fsolve(steady_state_func, init, xtol=1e-12, maxfev=2000)
        c_reactants[i] = sol[:lreactants]
        c_products[i] = sol[lreactants:-ladsorbed]
        theta[i] = sol[-ladsorbed:]
        fval[i] = steady_state_func(sol)
        init = sol

    return c_reactants, c_products, theta, fval, j
