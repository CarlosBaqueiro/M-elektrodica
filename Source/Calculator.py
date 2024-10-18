"""
M-elektrodica:
                Calculator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""
import numpy as np

from scipy.optimize import fsolve
from .Kpynetic import *

class Calculator:
    def __init__(self, data):

        self.c_reactants = np.zeros((len(data.E), len(data.reactants)))
        self.c_products = np.zeros((len(data.E), len(data.products)))
        self.theta = np.zeros((len(data.E), len(data.adsorbed)))
        self.fval = np.zeros((len(data.E), (len(data.reactants)+len(data.products)+len(data.adsorbed))))
        self.current = np.zeros(len(data.E))
        #Initialization
        c_reactants0 = np.ones(len(data.reactants))
        c_products0 = np.zeros(len(data.products))
        theta0 = np.zeros(len(data.adsorbed))
        initio = np.concatenate([c_reactants0, c_products0, theta0])

        for i in range(len(data.E)):
            solution = fsolve(self.steady_state, initio, args=(data.E[i], data))
            self.c_reactants[i], self.c_products[i], self.theta[i] = self.unzip_variables(solution, data)
            self.fval[i] = self.steady_state(solution, data.E[i], data)

    def unzip_variables(self, variables, data):
        c_reactants = variables[:len(data.reactants)]
        c_products = variables[len(data.reactants):-len(data.adsorbed)]
        theta = variables[-len(data.adsorbed):]
        return c_reactants, c_products, theta

    def steady_state(self, variables, V, data):
        c_reactants, c_products, theta = self.unzip_variables(variables, data)
        rhs = np.concatenate([(c_reactants - data.c0_reactants) * data.Fv / data.A,
                              (c_products - data.c0_products) * data.Fv / data.A,
                              np.zeros(len(theta))])

        return (dc_dt(aux(c_reactants, c_products, theta, V, data.Ga, data.DG, data.beta, data.ne,
                         data.nu, data.T), data.nux) - rhs)



