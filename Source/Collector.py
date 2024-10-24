import os
import numpy as np
import re

from .Tools import unit_conversion


class DataOperation:
    def __init__(self, operations_file):
        header, raw_data = Collector.raw_data(operations_file)
        self.variables_list = raw_data[:,header.index('Variables')]
        values = raw_data[:,header.index('Value')]
        # Kinetic model
        self.model = values[self.variables_list == 'Model']
        # Operations Conditions
            # Initialize cstr attribute
        cstr = values[self.variables_list == 'CSTR']
        if not cstr or cstr == 'False':
            self.cstr = False
        elif cstr == 'True':
            self.cstr = True
            self.Fv = float(values[self.variables_list == 'Volumetric flux'])

            # Temperature
        T = float(values[self.variables_list== 'Temperature'])
        uT = raw_data[:,header.index('Units')][self.variables_list== 'Temperature']
        self.T = unit_conversion('Temperature', T, uT, 'K')
        self.Ac = float(values[self.variables_list == 'Catalyst Active surface area'])

            #Potencial
        iE = float(values[self.variables_list == 'Initial potential'])
        fE = float(values[self.variables_list == 'Final potential'])
        hE = float(values[self.variables_list == 'Step potential'])
        self.potential = np.arange(iE, fE+hE, hE)

        #TODO: Add data recollection for more models and operations conditions

class DataSpecies:
    def __init__(self,species_file, operation):
        header, raw_data = Collector.raw_data(species_file)
        species_list = raw_data[:,header.index('Species')]
        index = header.index('RPACe')
        self.reactants = species_list[raw_data[:, index] == 'R'].tolist()
        self.products  = species_list[raw_data[:, index] == 'P'].tolist()
        self.adsorbed = species_list[raw_data[:, index] == 'A'].tolist()
        self.catalys   = species_list[raw_data[:, index] == 'C'].tolist()
        self.c0_reactants = np.array(raw_data[:, header.index('c0')][raw_data[:, index] == 'R'].astype(float))
        self.c0_products  = np.array(raw_data[:, header.index('c0')][raw_data[:, index] == 'P'].astype(float))
        self.list = self.reactants + self.products + self.adsorbed + self.catalys +['e-']
        if operation.model == 'Wang':
            self.DGads = np.array(raw_data[:,header.index('DGads')][raw_data[:, index] == 'A'].astype(float))

        #TODO: Add data recollection for more models

class DataReactions:
    def __init__(self, reaction_file, operation, species):
        header, raw_data = Collector.raw_data(reaction_file)
        self.list = raw_data[:, header.index('id')].tolist()
        self.beta = np.array(raw_data[:, header.index('Beta')].astype(float))
        self.nu = np.zeros((len(self.list), len(species.list)))

        for i in range(len(self.list)):
            r = raw_data[:, header.index('Reactions')][i]
            left, right = re.split(r'<->', r)

            speciesinreaction, stoichiometric = self.process_reaction(left, species.list)
            for specie, coeff in zip(speciesinreaction, stoichiometric):
                self.nu[i, species.list.index(specie)] = -coeff

            speciesinreaction, stoichiometric = self.process_reaction(right, species.list)
            for specie, coeff in zip(speciesinreaction, stoichiometric):
                self.nu[i, species.list.index(specie)] = coeff

        self.ne = self.nu[:,-1]
        self.nu = self.nu[:,:-1]
        self.nua = self.nu[:,len(species.reactants) + len(species.products): -len(species.catalys)]
        self.nux = self.nu[:,:len(species.reactants) + len(species.products) + len(species.adsorbed)]

        if operation.model == 'Wang':
            self.Ga = np.array(raw_data[:, header.index('Ga')].astype(float))
        #TODO: Add data recollection for more models

    @staticmethod
    def process_reaction(side, species_list):
        speciesinreaction =[]
        stoichiometric = []
        str = re.split(r'\s+', side)
        for s in str:
            s = s.strip()
            if s in ['', '+']:
                continue
            match = re.search(r"(\d+)?(.+)", s)
            if match:
                coefficient = int(match.group(1)) if match.group(1) else 1
                specie = match.group(2)
                if specie in species_list:
                    speciesinreaction.append(specie)
                    stoichiometric.append(float(coefficient))
                else:
                    raise ValueError(f"Process reaction error:  {specie} not found")
            else:
                raise ValueError(f"Process reaction error:  {str} not found")
        return speciesinreaction, stoichiometric


class Collector:
    def __init__(self, directory):
        self.operation = DataOperation(os.path.join(directory, 'operating.md'))
        self.species = DataSpecies(os.path.join(directory, 'species.md'), self.operation)
        self.reactions = DataReactions(os.path.join(directory, 'reactions.md'), self.operation, self.species)

    @staticmethod
    def raw_data(name_file):
        if not os.path.exists(name_file):
            raise FileNotFoundError(f"File was not found at the expected location: {name_file}")
        print(f'\t{name_file}')
        header = []
        raw_data = []
        with open(name_file, 'r') as f:
            lines = f.readlines()
            header = [col.strip() for col in lines[0].split('|')[1:-1]]
            for line in lines[2:]:
                if line.strip() and line[0] != '#':
                    entry_data = [col.strip() for col in line.split('|')[1:-1]]
                    raw_data.append(entry_data)
        raw_data = np.array(raw_data)
        return header, raw_data



