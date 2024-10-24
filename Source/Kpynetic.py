import numpy as np
from .Constants import *

#class Adsortion:
#    def __init__(self, data):
#        self.data = data
#        self.operation = data.operation
#        self.species = data.species
#        self.reactions = data.reactions
#    def Langmuir(self):
#        pass
#    def Frumkin(self):
#        pass

class Kpynetic:
    def __init__(self, data):
        self.data = data
        self.operation = data.operation
        self.species = data.species
        self.reactions = data.reactions

    def reaction_free_energy(self, upsilon, DG_formation):
        # The Gibbs free energy of the reaction for which this is written.
       return -(upsilon @ DG_formation)

    def reduction_reaction_free_energy(self, DGr, ne):
        return - ne * DGr

    def powerlaw(self, concentration, upsilon):
        return np.array([np.prod(concentration ** (-upsilon * (upsilon < 0)), axis=1),
                         -np.prod(concentration ** (upsilon * (upsilon > 0)), axis=1)])

    def rate(self, pre_exponential, k_rate, c_reactants, c_products, theta):
        concentrations = self.concentra(c_reactants, c_products, theta)
        rate = self.powerlaw(concentrations, self.reactions.nu)
        return pre_exponential*np.sum(k_rate*rate, axis=0)

    def dcdt(self, rate, upsilon: np.ndarray) -> np.ndarray:
        return np.dot(rate, upsilon)

    def emptysites(self, theta):
        return np.array([1-np.sum(theta[:-1])])
        #TODO: Modify for calculations with more than one catalyst

    def concentra(self, c_reactants, c_products, theta):
        return np.concatenate([c_reactants, c_products, theta, self.emptysites(theta)])

    def rate_constant (self, potential, ne, beta, T,
                       experimental_dependence,
                       chemical_dependence,
                       adsortion_dependence):
        return (experimental_dependence *
                np.exp (-(chemical_dependence + adsortion_dependence +
                          self.potential_dependence(potential, ne, beta)
                          ) / k_B / T))

    def experimental_dependence(self,  forward_constants , backward_constants):
        return np.array([forward_constants, backward_constants])

    def chemical_dependence(self, G_activation, DG_reaction, ne):
        DGrr = self.reduction_reaction_free_energy(DG_reaction, ne)
        return np.array([G_activation,  G_activation - DGrr])

    def potential_dependence(self, eta, ne, beta):
        return np.array([- ne * beta * eta,  ne * (1 - beta) * eta])

    def eta(self, potential, ne, DG_reaction):
        DG_reduction_reaction = self.reduction_reaction_free_energy(DG_reaction, ne)
        return abs(ne) * potential + DG_reduction_reaction

    def activation_free_energy(self, G_activation, DG_reaction, eta, ne, beta):
        DGrr = self.reduction_reaction_free_energy(DG_reaction, ne)
        return  G_activation - np.array([ne * beta * eta, DGrr - ne * (1 - beta) * eta])



class ExperimentalType(Kpynetic):
    def __init__(self, data, potential=0, c_reactants=None, c_products=None, theta=None):
        super().__init__(data)

        self.pre_exponential = 1
        self.eK = self.experimental_dependence(self.reactions.k_f, self.reactions.k_b)
        self.qK = 0 # Chemical dependence
        self.aK = 0 # Adsortion dependence

    def potential(self, potential, c_reactants, c_products, theta):
        self.k_rate = self.rate_constant(potential, self.reactions.ne, self.reactions.beta, self.operation.T,
                                         self.eK, self.qK, self.aK)
        self.v = self.rate(self.pre_exponential, self.k_rate, c_reactants, c_products, theta)

class WangType(Kpynetic):
    def __init__(self, data, potential=0, c_reactants=None, c_products=None, theta=None):
        super().__init__(data)

        self.pre_exponential = 1e3 / sum(self.reactions.ne)
        self.DG = self.reaction_free_energy(self.reactions.nua, self.species.DGads)
        self.eK = 1 # Experimental dependence
        self.qK = self.chemical_dependence(self.reactions.Ga, self.DG, self.reactions.ne)
        self.aK = 0 # Adsortion dependence
    def potential(self, potential, c_reactants, c_products, theta):
        self.k_rate = self.rate_constant(potential, self.reactions.ne, self.reactions.beta, self.operation.T,
                                         self.eK, self.qK, self.aK)
        self.v = self.rate(self.pre_exponential, self.k_rate, c_reactants, c_products, theta)





