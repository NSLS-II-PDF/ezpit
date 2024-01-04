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
    """
    database_atom_names = losa.load_atom_names(compton_aff_element_file)
    print('database_atom_names = ', database_atom_names)
    print all atom names in compton_element_only.txt
    database_atom_names =  ['H', 'He', 'Li', 'Be',,,,,,,,,,,'U']
    """
    with open(file_path, 'r') as f:
        atom_list = [line.strip() for line in f.readlines() if line.strip() != '']
    return atom_list


def load_scattering_factors(file_path):
    """
    database_scat_factors = losa.load_scattering_factors(compton_aff_parm_file)
    print('database_scat_factors = ', database_scat_factors)
    print all parameters for compton scattering form factor in compton_parameter_only.txt
    """
    with open(file_path, 'r') as f:
        data = []
        for line in f:
            data.append([float(val) for val in line.strip().split('\t')])
    return np.array(data)


def convert_atom_names(composition):
    """
    if composition = {'Co': 3, 'O': 4, 'P': 1},
    atom_names = losa.convert_atom_names(composition)
    print('atom_names = ', atom_names)
    --> atom_names = ['Co', 'Co', 'Co', 'O', 'O', 'O', 'O', 'P']
    Returns
    """
    atom_names = []
    for item in composition.items():
        key, value = item
        for j in range(value):
            atom_names.append(key)
    return atom_names


def get_scattering_factors(atom_names, database_atom_names,
                           database_scat_factors):
    """
    find Compton scattering parameters according to the given composition
    print('compton_scattering_factors ====== ', Compton_scattering_factors)
    # --> if composition is composition = {'Co':2, 'O':1, 'P':1}, get Compton scattering parameters for Co, O, P
    num_atom = len(atom_indices)
    # num_atom =  4, i.e, # N: total number of atoms in composition
    num_fact = len(compton_scattering_factors)
    # num_fact =  3, i.e. # how many different atoms in composition
    print('atomic_number = ', atomic_number)  # show atomic number of each atom in periodic table
    # --> if composition is composition = {'Co':2, 'O':1, 'P':1},  atomic_number =  [27, 8, 15]
    """
    scat_factors = []
    for atom in atom_names:
        if atom in database_atom_names:
            idx = database_atom_names.index(atom)
            scat_factors.append(database_scat_factors[idx])
        else:
            raise ValueError("There is no atom {0} in database".format(atom))
    return np.asarray(scat_factors)


#### added 12/18/2023######################################
def get_compton_scattering_factors(atom_names, database_atom_names,
                           database_scat_factors):
    """
    add compton scattering factor with give experimental q
    """

    scat_factors = []
    atomic_number = []
    for atom in atom_names:
        if atom in database_atom_names:
            idx = database_atom_names.index(atom)
            scat_factors.append(database_scat_factors[idx])
            atomic_number.append(idx+1)
        else:
            raise ValueError("There is no atom {0} in database".format(atom))
    return np.asarray(scat_factors), atomic_number
#### added 12/18/2023######################################


def group_atoms(atom_names):
    """
    Get unique atom names, their counts, and the index of atom in the unique
    name list. Results will be used by other functions.
            # if composition = {'Co':2, 'O':1, 'P':1},
            counter =  Counter({'Co': 2, 'O': 1, 'P': 1})
            atom_unique_names =  ['Co', 'O', 'P']
            atom_counts =  [2 1 1]  # "2" is for 'Co', "1" is for 'O', "1" is for 'P'
            atom_indices  =  [0 0 1 2] #
    """
    counter = Counter(atom_names)
    atom_uni_names = list(counter.keys())
    atom_counts = list(counter.values())
    atom_indices = [atom_uni_names.index(atom) for atom in atom_names]
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
