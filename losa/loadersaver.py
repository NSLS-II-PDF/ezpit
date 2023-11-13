import os
from collections import Counter
import numpy as np


def load_atom_name_positions(file_path):
    """
    Parameters
    ==========
    file_path: str
        Path to text file

    Returns
    =======

    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    atom_names = []
    atom_positions = []
    for line in lines:
        parts = line.split()
        atom_names.append(parts[0])
        atom_positions.append([float(x) for x in parts[1:]])
    atom_positions = np.array(atom_positions)
    return atom_names, atom_positions


def load_atom_names(file_path):
    with open(file_path, 'r') as f:
        atom_list = [line.strip() for line in f.readlines() if line.strip() != '']
    return atom_list

def load_scattering_factors(file_path):
    with open(file_path, 'r') as f:
        data = []
        for line in f:
            data.append([float(val) for val in line.strip().split('\t')])
    return np.array(data)

def convert_atom_names(composition):
    atom_names = []
    for item in composition.items():
        key, value = item
        for j in range(value):
            atom_names.append(key)
    return atom_names

def get_scattering_factors(atom_names, database_atom_names,
                           database_scat_factors):
    scat_factors = []
    for atom in atom_names:
        if atom in database_atom_names:
            idx = database_atom_names.index(atom)
            scat_factors.append(database_scat_factors[idx])
        else:
            raise ValueError("There is no atom {0} in database".format(atom))
    return np.asarray(scat_factors)


def group_atoms(atom_names):
    """
    Get unique atom names, their counts, and the index of atom in the unique
    name list. Results will be used by other functions.
    """
    counter = Counter(atom_names)
    #print('conter = ', counter)
    atom_uni_names = list(counter.keys())
    #print('atom_uni_names = ', atom_uni_names)
    atom_counts = list(counter.values())
    #print('atom_counts = ', atom_counts)
    atom_indices = [atom_uni_names.index(atom) for atom in atom_names]
    #print('atom_indices = ', atom_indices)
    return atom_uni_names, np.asarray(atom_counts), np.asarray(atom_indices)


def make_folder(file_path):
    """
    Create a folder for saving file if the folder does not exist. This is a
    supplementary function for savers.
    Parameters
    ----------
    file_path : str
        Path to a file.
    """
    file_base = os.path.dirname(file_path)
    if not os.path.exists(file_base):
        try:
            os.makedirs(file_base, exist_ok=True)
        except OSError:
            raise ValueError("Can't create the folder: {}".format(file_base))


def save_txt(filename, q_Iq):
    make_folder(filename)
    np.savetxt(filename, q_Iq)
