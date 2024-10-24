"""
M-elektrodica:
                Calculator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""
# Importar las bibliotecas necesarias
import numpy as np
from scipy.optimize import fsolve
from .Kpynetic import WangType
from .Grapher import Grapher

class BaseConcentration:
    def __init__(self, data, Kpy):
        self.data = data
        self.operation = data.operation
        self.species = data.species
        self.reactions = data.reactions
        self.Kpy = Kpy

    def simulator(self, Kpy):
        self.c_reactants = np.zeros((len(self.operation.potential), len(self.species.reactants)))
        self.c_products = np.zeros((len(self.operation.potential), len(self.species.products)))
        self.theta = np.zeros((len(self.operation.potential), len(self.species.adsorbed)))
        self.j = np.zeros(len(self.operation.potential))
        self.fval, initio = self.inicialization(self.operation, self.species)
        for i in range(len(self.operation.potential)):
            potential = self.operation.potential[i]
            solution = fsolve(self.steady_state, initio,
                                  args=(potential, Kpy, self.data), xtol=1e-12, maxfev=2000)
            self.c_reactants[i], self.c_products[i], self.theta[i] = self.unzip_variables(solution, self.species)
            self.fval[i] = self.steady_state(solution, potential, Kpy, self.data)
            self.j[i] = self.current(solution, potential, Kpy, self.data)
            initio = solution
        return self

    def inicialization(self, operation, species):
        raise NotImplementedError("The method inicialization must be implemented by the subclass")

    def unzip_variables(self, variables, species):
        raise NotImplementedError("The method unzip_variables must be implemented by the subclass")

    def ride_hand_side(self, c_reactants, c_products, theta, data):
        raise NotImplementedError("The method unzip_variables must be implemented by the subclass")

    def steady_state(self, variables, potential, Kpy, data):
        c_reactants, c_products, theta = self.unzip_variables(variables, data.species)
        rhs = self.ride_hand_side(c_reactants, c_products, theta, data)
        Kpy.potential(potential, c_reactants, c_products, theta)
        return Kpy.dcdt(Kpy.v, data.reactions.nux) - rhs

    def current(self, variables, potential, Kpy, data):
        c_reactants, c_products, theta = self.unzip_variables(variables, data.species)
        Kpy.potential(potential, c_reactants, c_products, theta)
        return np.dot(data.reactions.ne, Kpy.v)


class StaticConcentration(BaseConcentration):
    def inicialization(self, operation, species):
        fval = np.zeros((len(operation.potential), len(species.adsorbed)))
        theta0 = np.zeros(len(species.adsorbed))
        initio = np.concatenate([theta0])
        return fval, initio

    def unzip_variables(self, variables, species):
        c_reactants = species.c0_reactants
        c_products = species.c0_products
        theta = variables
        return c_reactants, c_products, theta

    def ride_hand_side(self, c_reactants, c_products, theta, data):
        return np.zeros(len(theta))

class DynamicConcentration(BaseConcentration):
    def inicialization(self, operation, species):
        fval = np.zeros((len(operation.potential),
                         len(species.reactants) + len(species.products) + len(species.adsorbed)))
        c_reactants0 = np.ones(len(species.reactants))
        c_products0 = np.zeros(len(species.products))
        theta0 = np.zeros(len(species.adsorbed))
        initio = np.concatenate([c_reactants0, c_products0, theta0])
        return fval, initio

    def unzip_variables(self, variables, species):
        c_reactants = variables[:len(species.reactants)]
        c_products = variables[len(species.reactants):-len(species.adsorbed)]
        theta = variables[-len(species.adsorbed):]
        return c_reactants, c_products, theta

    def ride_hand_side(self, c_reactants, c_products, theta, data):
        return np.concatenate([(c_reactants - data.species.c0_reactants) * data.operation.Fv / data.operation.Ac,
                               (c_products - data.species.c0_products) * data.operation.Fv / data.operation.Ac,
                               np.zeros(len(theta))])

class Calculator:
    def __init__(self, data):
        self.data = data
        self.operation = data.operation

        if self.operation.model == 'Wang':
            self.Kpy = WangType(self.data)

        if self.operation.cstr:
            solver = DynamicConcentration(self.data, self.Kpy)
        else:
            solver = StaticConcentration(self.data, self.Kpy)

        self.results = solver.simulator(self.Kpy)
        Grapher(self.data, self.results)

