import os
import numpy as np
import re

class RecuperaData:
    def __init__(self, directory):
        species_file = os.path.join(directory, 'species_list.md')
        reaction_file = os.path.join(directory, 'reactions_list.md')
        operations_file = os.path.join(directory, 'operation_conditions.md')

        self.reactants, self.c0_reactants, self.products, self.c0_products, self.adsorbed = self.species(species_file)
        self.Ga, self.DG, self.beta, self.nu = self.reactions(reaction_file, self.reactants, self.products,
                                                              self.adsorbed)
        self.ne = self.nu[:, -1]
        self.nux = self.nu[:, :(len(self.reactants) + len(self.products) + len(self.adsorbed))]
        self.nu = self.nu[:, :(len(self.reactants) + len(self.products) + len(self.adsorbed)) + 1]

        variables_list = ['T', 'Fv', 'A', 'E_i', 'E_f', 'E_h']
        variables = self.op_variables(operations_file, variables_list)
        self.T = float(variables['T'])
        self.Fv = float(variables['Fv'])
        self.A = float(variables['A'])
        self.E = np.arange(float(variables['E_i']), float(variables['E_f'])+float(variables['E_h']), float(variables['E_h']))

    @staticmethod
    def species(file_species):
        reactants = []
        c0_reactants = []
        products = []
        c0_products = []
        adsorbed = []
        with open(file_species, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                if line.strip() and line[0] != '#':
                    columns = [col.strip() for col in line.split('|')[1:-1]]
                    if columns[1] == 'R':
                        reactants.append(columns[0])
                        c0_reactants.append(float(columns[2]))
                    elif columns[1] == 'P':
                        products.append(columns[0])
                        c0_products.append(float(columns[2]))
                    elif columns[1] == 'A':
                        adsorbed.append(columns[0])
        c0_reactants = np.array(c0_reactants)
        c0_products = np.array(c0_products)
        return reactants, c0_reactants, products, c0_products, adsorbed

    @staticmethod
    def process_species(side, species_list):
        u = np.zeros(len(species_list))
        species = re.split(r'\s+', side)
        for s in species:
            s = s.strip()
            if s in ['', '+']:
                continue
            match = re.search(r"(\d+)?(.+)", s)
            if match:
                coefficient = int(match.group(1)) if match.group(1) else 1
                compound = match.group(2)
                if compound in species_list:
                    u[species_list.index(compound)] = coefficient
                else:
                    raise ValueError(f" {compound} error")
            else:
                raise ValueError(f"No se pudo hacer match con el compuesto {s}")
        return u

    @staticmethod
    def reactions(file_reactions, reactants, products, adsorbed):
        species_list = reactants + products + adsorbed + ['Pt', 'H+', 'e-']
        Ga = []
        DG = []
        beta = []
        nu = []
        with open(file_reactions, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                if line.strip() and line[0] != '#':
                    columns = [col.strip() for col in line.split('|')[1:-1]]
                    left, right = re.split(r'<->', columns[0])
                    R = RecuperaData.process_species(left, species_list)
                    P = RecuperaData.process_species(right, species_list)
                    nu.append(P - R)
                    Ga.append(float(columns[1]))
                    DG.append(float(columns[2]))
                    beta.append(float(columns[3]))
        Ga = np.array(Ga)
        DG = np.array(DG)
        beta = np.array(beta)
        nu = np.array(nu)
        return Ga, DG, beta, nu

    @staticmethod
    def op_variables(file_operations, variables_list):
        variables = {}
        with open(file_operations, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                if line.strip() and line[0] != '#':
                    columns = [col.strip() for col in line.split('|')[1:-1]]
                    variable = columns[1]
                    valor = columns[2]
                    if variable in variables_list:
                        variables[variable] = valor
        return variables
