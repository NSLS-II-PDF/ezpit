import matplotlib.pyplot as plt
import losa.loadersaver as losa
import proc.processing as proc
import timeit
import numpy as np
import sys

input_base = "C:/Users/gkwon/Pycharmprojects/ezpit/data/"
compton_aff_element_file = input_base + 'compton_element_only.txt'
compton_aff_parm_file = input_base + 'compton_parameter_only.txt'
compton_atomnumber_file = input_base + 'compton_atomicnumber.txt'

#parameters
wavelength = 0.1665
alpha = 3  #2 or 3 can be used.
composition = {'Co':2, 'O':2, 'P':1}
qmin = 0.6
qmax = 24
qstep = 0.01

atom_names = losa.convert_atom_names(composition)
# Database
# #load all of atom's name and its parameters for compton scattering form factor
compton_atom_names = losa.load_atom_names(compton_aff_element_file)
compton_scat_parms = losa.load_scattering_factors(compton_aff_parm_file)
atom_unique_names, atom_counts, atom_indices = losa.group_atoms(atom_names)
print('atom_unique_names = ',atom_unique_names)
print('atom_counts = ',atom_counts)
print('atom_indices = ',atom_indices)

#laod only parameters of compton scattering form factor for the given composition:
#both (get_scattering_factors or get_compton_scattering_factors) are working
#compton_scat_form_factor, atomic_number = losa.get_scattering_factors(atom_unique_names,
#                                                           compton_atom_names,compton_scat_parms) #[1]

compton_scat_form_factor, atomic_number = losa.get_compton_scattering_factors(atom_unique_names,
                                                                    compton_atom_names,compton_scat_parms) #[1]

list_q, list_compton_scat = proc.compton_calc_exp(atom_indices, compton_scat_parms, compton_scat_form_factor,
                            atomic_number, qmin=qmin, qmax=qmax, qstep=qstep, wavelength=wavelength, alpha=alpha)

np.savetxt(input_base + 'list_compton_scat.chi', np.column_stack(([list_q, list_compton_scat]))) # or use "list(zip(r, Gr)))"

plt.figure(0)
plt.plot(list_q, list_compton_scat, label='Compton_scat_pattern')
plt.grid()
plt.legend()
plt.show()
