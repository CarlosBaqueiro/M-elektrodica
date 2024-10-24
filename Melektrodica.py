"""
M-elektrodica:
                Simulator
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""

import os
from Source import *
begin()
#pakage = 'Sheng2020'
pakage = 'Wang2007Hydrogen'

if __name__ == '__main__':
    #ValidationTest()
    run = True
    test = False
    if run:
        directory = os.path.join('Examples', pakage)
        if not os.path.exists(directory):
            print(f"The directory {directory} does not exist.")
        else:
            print(f"Directory: {directory}")
            data = Collector(directory)
            results = Calculator(data)
            if test:
                Tester(data, results)

