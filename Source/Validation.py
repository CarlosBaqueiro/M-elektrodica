import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#sys.exit()


from .Kpynetic import *
from .Collector import *
from .Calculator import *
from. Tools import showme


class Hydrogen:
    def __init__(self):
        directory = os.path.join('Examples', 'Wang2007Hydrogen')
        data = Collector(directory)
        self.operation = data.operation
        self.species = data.species
        self.reactions = data.reactions

        # Expected data
        self.name_reactions = ['DA', 'OA', 'OD']
        self.nu = np.array([[-1.0, 0.0, 2.0, -2.0], [-1.0, 1.0, 1.0, -1.0], [ 0.0, 1.0, -1.0, 1.0]])
        self.nua = [[2.0],[ 1.0],[-1.0]]
        self.nui = self.nu[:,-2:]
        self.ne = np.array ([0.0, 1.0, 1.0])
        self.Ga = np.array([196, 294, 48]) * 1e-3
        self.DGads = np.array([75]) * 1e-3
        self.beta = np.array ([0.0, 0.5, 0.5])

        ValidationTest.data_validation(self.nu, self.reactions.nu, 'Coefficient matrix')
        ValidationTest.data_validation(self.ne, self.reactions.ne, 'Electrons transferred vector')
        ValidationTest.data_validation(self.Ga, self.reactions.Ga, 'Activation energy vector')
        ValidationTest.data_validation(self.DGads, self.species.DGads, 'Free energy vector')


        self.results = Calculator(data)
        data.operation.cstr = True
        data.operation.Fv = 1
        data.operation.Ac = 1
        #data.species.c0_reactants = 0
        self.results_cstr = Calculator(data)
        # Expected results
        # Wang 2007 Hydrogen, Figure 2
        theta_exp =np.array([[-0.00609637,  0.06183049],
                               [ 0.00660171,  0.03845448],
                               [ 0.01890746,  0.02559091],
                               [ 0.03082574,  0.01693716],
                               [ 0.04312665,  0.01037622],
                               [ 0.05542675,  0.00486571],
                               [ 0.06772363,  0.00355696],
                               [ 0.07963867,  0.00089504],
                               [ 0.09193232,  0.00199795],
                               [ 0.1042292 ,  0.00068919],
                               [ 0.11652526,  0.00043087],
                               [ 0.12882133,  0.00017255],
                               [ 0.14111659,  0.00096467],
                               [ 0.15341266,  0.00070635],
                               [ 0.16532448,  0.0004561 ],
                               [ 0.1780048 ,  0.0001897 ],
                               [ 0.19030006,  0.00098182],
                               [ 0.20144337,  0.00074772]])

        self.potential = theta_exp[:,0]
        self.theta = theta_exp[:,1]

        # Wang 2007 Hydrogen, Equation 33
        self.theta_eq = np.zeros(len(self.potential))
        self.DG = Kpynetic.reaction_free_energy(self.nua, self.DGads)
        #DG = self.nua @ self.DGads #ValidationTest.reaction_free_energy(self.nua, self.DGads)
        #print(DG, '\n')

        #diff = np.zeros((len(self.potential), 2 * len(self.name_reactions)))
        for i in range(len(self.potential)):
            E = self.potential[i]
            self.DGa = Kpynetic.activation_free_energy(self.Ga, self.DG, E, self.ne, self.beta)
            self.g = np.exp(-self.DGa/ k_B / self.operation.T)
            #g = ValidationTest.gibbs_free(self.Ga, self.DGads, E, self.nua, self.ne, self.beta, self.operation.T)
            #diff[i] = (self.g - g).flatten()
            g=self.g
            # Wang 2007 Hydrogen, Equation 34-37
            A = 2*g[0, 0]-2*g[1, 0]
            B = -4*g[0, 0] - g[0, 1] - g[1, 1] - g[0, 2] - g[1, 2]
            C = 2*g[0, 0] + g[1, 1] + g[1, 2]
            if A == 0:
                self.theta_eq[i] = -C/B
            else:
                self.theta_eq[i] = (-B - np.sqrt(B**2 - 4*A*C))/(2*A)



        #for h in range(len(self.name_reactions)):
        #    plt.plot(self.potential, diff[:,h], label=f'{self.name_reactions[h]}_f')
        #    plt.plot(self.potential, diff[:,h+len(self.name_reactions)], label=f'{self.name_reactions[h]}_b')
        #    plt.legend( )
        #    plt.show()
        ValidationTest.adsorbed_validation(self.potential, self.theta, self.theta_eq,
                                           self.results, self.operation, self.results_cstr)
        #ValidationTest.current_validation(self.potential, self.theta, self.theta_eq,
        #                                   self.results, self.operation, self.results_cstr)
class ValidationTest:
    def __init__(self):
        Hydrogen()

    @staticmethod
    def gibbs_free(Ga, DGads, eta, nua, ne, beta, T):
        DG = Ga - np.array([ne * beta * eta, nua@DGads - ne * (1 - beta) * eta])
        g = np.exp(-DG / k_B / T)
        return g

    @staticmethod
    def data_validation(expected, collect, variable):
        try:
            if not np.array_equal(expected, collect):
                raise ValueError(f"Validation Error: {variable}")
            print(f"Validation Success: {variable}.")
        except ValueError as ve:
            print(ve)
            raise
    @staticmethod
    def adsorbed_validation(potential, theta, theta_eq, results, operation, results_cstr):
        plt.plot(potential, theta_eq, label=r'Wang2007 (Equations 34-37)')
        plt.plot(potential, theta, marker='o',linestyle='none', label=r'Wang2007 (Fig. 2)')
        plt.plot(operation.potential, results.theta, label=r'Results')
        plt.plot(operation.potential, results_cstr.theta, label=r'Results CSTR')
        plt.yscale('log')
        plt.xlabel(r'$\eta$ [V]')
        plt.ylabel(r'$\theta_{H*}$')
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()
        plt.legend( )
        plt.show()

