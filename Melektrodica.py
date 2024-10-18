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
from Source.Tests import *

#pakage = 'Sheng2020'
pakage = 'Wang2007Hydrogen'

if __name__ == '__main__':
    directory = os.path.join('Examples', pakage)
    test = True
    if not os.path.exists(directory):
        print(f"The directory {directory} does not exist.")
    else:
        data = RecuperaData(directory)
        results = Calculator(data)
        if test:
            Test(data, results)
        else:
            Grapher(data, results)