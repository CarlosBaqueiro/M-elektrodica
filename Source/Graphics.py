"""
M-elektrodica:
                Graphics functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""


import matplotlib.pyplot as plt

def fval(E, f, variables):
    plt.plot(E, f, linestyle='-.', label=variables)
    plt.xlabel('Potential [V]')
    plt.ylabel(r'$\frac{\partial \theta_i}{\partial t}$')
    ax = plt.gca()
    ax.yaxis.label.set_size(15)
    plt.yscale('log')
    plt.grid(visible=True, which='both', axis='both',color='grey', linestyle='-', linewidth='0.2')
    plt.minorticks_on()
    plt.legend()
    plt.tight_layout()
    plt.show()