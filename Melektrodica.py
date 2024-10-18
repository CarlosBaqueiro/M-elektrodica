"""
M-elektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import os
from Source.Calculator import *


if __name__ == '__main__':
    directory = os.path.join('Examples', 'Sheng2020')

    if not os.path.exists(directory):
        print(f"El directorio {directory} no existe.")
    else:
        Calculator(directory)




