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
        operation = data.operation
        species = data.species
        reactions = data.reactions
        self.Kpy = Kpynetic(operation, species, reactions)
        self.c_reactants = np.zeros((len(operation.potential), len(species.reactants)))
        self.c_products = np.zeros((len(operation.potential), len(species.products)))
        self.theta = np.zeros((len(operation.potential), len(species.adsorbed)))


        self.j = np.zeros(len(operation.potential))
        #Initialization

        theta0 = np.zeros(len(species.adsorbed))
        catalys0 = np.ones(len(species.catalys))


        if operation.cstr == False:
            self.fval = np.zeros((len(operation.potential), len(species.adsorbed)))
            initio = np.concatenate([theta0])
            for i in range(len(operation.potential)):
                solution = fsolve(self.steady_state_const, initio, args=(operation.potential[i], species, reactions, operation),
                                  xtol=1e-20, maxfev=2000)
                self.theta[i] = solution
                self.c_reactants[i] = species.c0_reactants
                self.c_products[i] = species.c0_products
                self.fval[i] = self.steady_state_const(solution, operation.potential[i], species, reactions, operation)
                self.j[i] = self.current_const(solution, operation.potential[i], species, reactions, operation)
                initio = solution
        else:
            c_reactants0 = np.ones(len(species.reactants))
            c_products0 = np.zeros(len(species.products))
            self.fval = np.zeros((len(operation.potential), len(species.list)-2))
            initio = np.concatenate([c_reactants0, c_products0, theta0])
            for i in range(len(operation.potential)):
                solution = fsolve(self.steady_state, initio, args=(operation.potential[i], species, reactions, operation),
                                  xtol=1e-20, maxfev=2000)
                self.c_reactants[i], self.c_products[i], self.theta[i] = self.unzip_variables(solution, species)
                self.fval[i] = self.steady_state(solution, operation.potential[i], species, reactions, operation)
                self.j[i] = self.current(solution, operation.potential[i], species, reactions, operation)
                initio = solution

    def aux_const(self, theta, V, species, reactions, operation):
        c = self.Kpy.concentra(species.c0_reactants, species.c0_products, theta)
        w = self.Kpy.wang(V, species, reactions, operation)
        self.rate = Kpynetic.rate(w.preexponential, w.rate_constants, self.Kpy.powerlaw(c, reactions.nu))
        return self.rate

    def steady_state_const(self, theta, V, species, reactions, operation):
        return self.Kpy.dcdt(self.aux_const(theta, V, species, reactions, operation), reactions.nua)

    def current_const(self, theta, V, species, reactions, operation):
        return np.dot(reactions.ne, self.aux_const(theta, V, species, reactions, operation))

    def unzip_variables(self, variables, species):
        c_reactants = variables[:len(species.reactants)]
        c_products = variables[len(species.reactants):-len(species.adsorbed)]
        theta = variables[-len(species.adsorbed):]
        return c_reactants, c_products, theta

    def steady_state(self, variables, V, species, reactions, operation):
        self.aux(variables, V, species, reactions, operation)
        return self.Kpy.dcdt(self.rate, reactions.nux) - self.rhs

    def aux(self, variables, V, species, reactions, operation):
        c_reactants, c_products, theta = self.unzip_variables(variables, species)
        c = self.Kpy.concentra(c_reactants, c_products, theta)
        w = self.Kpy.wang(V, species, reactions, operation)
        self.rate = Kpynetic.rate(w.preexponential, w.rate_constants, self.Kpy.powerlaw(c, reactions.nu))
        self.rhs = np.concatenate([(c_reactants - species.c0_reactants) * operation.Fv / operation.Ac,
                                   (c_products - species.c0_products) * operation.Fv / operation.Ac,
                                   np.zeros(len(theta))])
        return

    def current(self, variables, V, species, reactions, operation):
        self.aux(variables, V, species, reactions, operation)
        return np.dot(reactions.ne, self.rate) * F * k_B * operation.T / h
