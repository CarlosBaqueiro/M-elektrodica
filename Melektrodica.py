"""
M-elektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import os
from Source.Calculator import *
from Source.Data_collection import *
from Source.Graphics import *
from Source.Tests import *
from Source.Validation_test import *
#pakage = 'Sheng2020'
pakage = 'Wang2007Hydrogen'

if __name__ == '__main__':
    ValidationTest()
    run = True
    test = False
    if run:
        directory = os.path.join('Examples', pakage)
        if not os.path.exists(directory):
            print(f"The directory {directory} does not exist.")
        else:
            data = Collector(directory)
            results = Calculator(data)
            Grapher(data, results)
            if test:
                Tester(data, results)

