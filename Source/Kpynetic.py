import numpy as np
from .Constants import *

class Kpynetic:
    def __init__(self, operations, species, reactions):
        self.operations = operations
        self.species = species
        self.reactions = reactions

    def experimental(self, eta, reactions):
        self.rate_constants = np.array([reactions.k_f, reactions.k_b])

    def wang(self, eta, species, reactions, operation):
        self.DG  = Kpynetic.reaction_free_energy(reactions.nua, species.DGads)
        self.DGa = Kpynetic.activation_free_energy(reactions.Ga, self.DG, eta, reactions.ne, reactions.beta)
        self.rate_constants = np.exp(-self.DGa/ k_B / operation.T)
        self.preexponential = 1e3/2
        return self

    @staticmethod
    def reaction_free_energy(upsilon, DG_formation):
        # The Gibbs free energy of the reaction for which this is written.
        return -(upsilon @ DG_formation)

    @staticmethod
    def reduction_reaction_free_energy(DGr, ne):
        return - ne * DGr

    @staticmethod
    def activation_free_energy(G_activation, DG_reaction, eta, ne, beta):
        DGrr = Kpynetic.reduction_reaction_free_energy(DG_reaction, ne)
        return  G_activation - np.array([ne * beta * eta, DGrr - ne * (1 - beta) * eta])

    @staticmethod
    def rate_constant (DG_activation, T):
        return np.exp(-DG_activation / k_B / T)

    @staticmethod
    def powerlaw(concentration, upsilon):
        return np.array([np.prod(concentration ** (-upsilon * (upsilon < 0)), axis=1),
                         -np.prod(concentration ** (upsilon * (upsilon > 0)), axis=1)])
    @staticmethod
    def rate (pre_exponential, k_constants, powerlaw) -> np.ndarray:
        return pre_exponential*np.sum(k_constants*powerlaw, axis=0)

    @staticmethod
    def dcdt (rate, upsilon: np.ndarray) -> np.ndarray:
        return np.dot(rate, upsilon)

    @staticmethod
    def emptysites (theta, upsilon):
        #return np.array[]
        pass

    @staticmethod
    def concentra (c_reactants, c_products, theta):
        #TODO: Modify for calculations with more than one catalyst
        return np.concatenate([c_reactants, c_products, theta, np.array([1-np.sum(theta[:-1])])])
