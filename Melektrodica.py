"""
M-elektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

#from Source.include import *
import os
from Source import Recupera_data, steady_state_function



def calculator (directory):
    data = Recupera_data.keyData(directory)
    result = steady_state_function.SteadyState(data)

    return result

if __name__ == '__main__':
    directory = os.path.join('Examples', 'Sheng2020')

    if not os.path.exists(directory):
        print(f"El directorio {directory} no existe.")
    else:
        print(f"El directorio {directory} existe.")
        print(calculator(directory))


