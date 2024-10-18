"""
M-elektrodica:
                Graphics functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""


import matplotlib.pyplot as plt
import numpy as np


class Grapher:
    def __init__(self, data, results):

        self.fval(data.E, results.fval, data.reactants + data.products + data.adsorbed)
        self.current(data.E, results.j )

    def graph_results(self, data, results):
        self.species(data, results, 'inSolution', login=False)
        self.species(data, results, 'coverages', login=True)

    def fval(self, E, f, variables):
        plt.plot(E, f, linestyle='-.', label=variables)
        plt.xlabel('Potential [V]')
        plt.ylabel(r'$\frac{\partial \theta_i}{\partial t}$')
        ax = plt.gca()
        ax.yaxis.label.set_size(15)
        # plt.yscale('log')
        plt.grid(True, which='both', axis='both', color='grey', linestyle='-', linewidth='0.2')
        plt.minorticks_on()
        plt.legend()
        plt.tight_layout()
        plt.show()

    def species(self, data, results, label, login=False):
        if label == 'coverages':
            c_species = results.theta
            legends = data.adsorbed
            plt.ylabel(r'$\theta_i$')
        elif label == 'inSolution':
            c_species = np.concatenate([results.c_reactants, results.c_products], axis=1)
            legends = data.reactants + data.products
            plt.ylabel(r'$c_i\ mol/L$')

        plt.plot(data.E, c_species, label=legends)
        plt.xlabel('Potential [V]')
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.minorticks_on()
        plt.legend(loc='lower right')
        if login: plt.yscale('log')
        plt.tight_layout()
        plt.show()

    def current(self, E, j):
        plt.plot(j, E)
        plt.xlabel('Current density [A/cm2]')
        plt.ylabel('Ptencial [V]')
        #plt.xscale('log')
        plt.tight_layout()
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.show()