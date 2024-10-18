"""
M-elektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import os
from Source.Calculator import *
from Source.Recupera_data import *
from Source.Graphics import *


if __name__ == '__main__':
    directory = os.path.join('Examples', 'Sheng2020')

    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
    else:
        data = RecuperaData(directory)
        results = Calculator(data)
        Grapher(data, results)



