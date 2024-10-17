"""
M-elektrodica:
                Tools functions
@author : Carlos Baqueiro Basto

Created on Julio 2024
"""
import numpy as np
import re

def species_data(file_species):
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
                raise ValueError(f"El compuesto {compound} no está en la lista de especies")
        else:
            raise ValueError(f"No se pudo hacer match con el compuesto {s}")
    return u


def reaction_data(file_reactions, reactants, products, adsorbed):
    species_list = reactants + products + adsorbed + ['Pt', 'H+', 'e-']
    Ga = []
    DG = []
    nu = []
    with open(file_reactions, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:
            if line.strip() and line[0] != '#':
                columns = [col.strip() for col in line.split('|')[1:-1]]
                Ga.append(columns[0])
                DG.append(columns[1])
                left, right = columns[2].split(r'<-->')
                R = process_species(left, species_list)
                P = process_species(right, species_list)
                nu.append(P - R)
    Ga = np.array(Ga)
    DG = np.array(DG)
    nu = np.array(nu)
    return Ga, DG, nu