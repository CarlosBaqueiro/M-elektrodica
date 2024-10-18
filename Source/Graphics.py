"""
M-elektrodica:
                Graphics functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""


import matplotlib.pyplot as plt
class Grapher:
    def __init__(self, data, results):

        self.fval(data.E, results.fval, data.reactants + data.products + data.adsorbed)
        self.coverages(data.E, results.theta, data.adsorbed)
        self.current(results.current, data.E)

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

    def coverages(self, E, theta, adsorbed):
        plt.plot(E, theta, label=adsorbed)
        #plt.plot(E, (1-np.sum(theta*sites, axis=1)), label=r'$1-\sum \theta_i$')
        #plt.yscale('log')
        plt.xlabel('Potential [V]')
        plt.ylabel(r'$\theta_i$')
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.minorticks_on()
        plt.tight_layout()
        plt.legend(loc='lower right')
        plt.show()

    def current(self, E, j):
        plt.plot(j, E)
        plt.xlabel('Current density [A/cm2]')
        plt.ylabel('Ptencial [V]')
        #plt.xscale('log')
        plt.tight_layout()
        plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-',linewidth='0.2')
        plt.show()