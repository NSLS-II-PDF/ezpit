import matplotlib.pyplot as plt
import losa.loadersaver as losa
import proc.processing as proc
import timeit
import numpy as np
import sys

#Scattering form factor data file
input_base = "C:/Users/gkwon/Pycharmprojects/ezpit/data/"
aff_element_file = input_base + 'aff_elementonly.txt'
aff_parm_file = input_base + 'aff_parmonly.txt'

#Compton scattering data file
input_base = "C:/Users/gkwon/Pycharmprojects/ezpit/data/"
compton_aff_element_file = input_base + 'compton_element_only.txt'
compton_aff_parm_file = input_base + 'compton_parameter_only.txt'
compton_atomnumber_file = input_base + 'compton_atomicnumber.txt'

#load experimental Iq (expqIq_data and background I(q) (bkgqIq_data)
#The skipping header for 'xxx.chi' file is located within the "calculate_expSq()' function".
#expqIq_data = 'C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5.chi'
#bkgqIq_data = 'C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_0p7cap_Nsum6.chi'

expIq_input_base = "C:/Users/gkwon/PycharmProjects/ezpit/dev/Test_qmax_file/Chi_1file/"
bkgqIq_input_base = "C:/Users/gkwon/PycharmProjects/ezpit/dev/Test_qmax_file/background_Chifile_1file/"
expqIq_data = expIq_input_base + 'A_CoPiITOglass_02142024_1-test_20240219-154304_696917_primary-dk_sub_image-0.chi'
bkgqIq_data = bkgqIq_input_base + 'A_emptyquartzcap_0p5_20240219-122602_4e50c5_primary-dk_sub_image-0.chi'

#Parameters
composition = {'Co':1, 'O':1, 'P': 1}
qmin = 0.6
qmax = 24
qstep = 0.01
background_scale = 0.013
qdamp = 0
poly_order = 11
rmin = 0
rmax = 20
rstep = 0.01
wavelength = 0.1665
alpha = 3

################################### Compton scattering parts  added 12/18/2023#####################
atom_names = losa.convert_atom_names(composition)
# Database
compton_atom_names = losa.load_atom_names(compton_aff_element_file)
#load all of parameters for compton scattering form factor
compton_scat_parms = losa.load_scattering_factors(compton_aff_parm_file)
atom_unique_names, atom_counts, atom_indices = losa.group_atoms(atom_names)

#laod only parameters of compton scattering form facotr for the given composition
compton_scat_form_factor, atomic_number = losa.get_compton_scattering_factors(atom_unique_names,
                                                                      compton_atom_names,compton_scat_parms) #[1]

list_q, list_compton_scat = proc.compton_calc_exp(atom_indices, compton_scat_parms, compton_scat_form_factor,
                                atomic_number, qmin=qmin, qmax=qmax, qstep=qstep, wavelength=wavelength, alpha=alpha)

np.savetxt(input_base + 'list_compton_scat.chi', np.column_stack(([list_q, list_compton_scat]))) # or use "list(zip(r, Gr)))"

plt.figure(10)
plt.plot(list_q, list_compton_scat, label='compton_scat')
plt.xlabel('q (1/A)')
plt.ylabel('Compton scatt(q)')
plt.grid()
plt.legend()
################################### Compton scattering parts  added 12/18/2023#####################


# Database
database_atom_names = losa.load_atom_names(aff_element_file)
database_scat_factors = losa.load_scattering_factors(aff_parm_file)

atom_unique_names, atom_counts, atom_indices = losa.group_atoms(atom_names)

scattering_factors = losa.get_scattering_factors(atom_unique_names, database_atom_names,
                                                database_scat_factors)
q, Iq, scaled_expIq, list_scaled_bkgIq, list_Sq, Sq, Fq, mean_sq_fi, sq_mean_fi = proc.calculate_expSq(atom_indices,
            scattering_factors, expqIq_data, bkgqIq_data, qmin=qmin, qmax=qmax, qstep=qstep,
            background_scale=background_scale, poly_order=11, return_Iq=False)

#load test results in xPDFsuite with the same data plotted here.
test_xpdfsuite_qsq = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-Nohdeader.sq", skiprows=26)
test_xpdfsuite_qfq = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-NOheader.fq", skiprows=26)
test_xpdfsuite_rgr = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-Noheader.gr", skiprows=26)

np.savetxt(input_base + 'square_mean_fi.chi', np.column_stack(([list_q, sq_mean_fi]))) # or use "list(zip(r, Gr)))"
np.savetxt(input_base + 'mean_square_fi.chi', np.column_stack(([list_q, mean_sq_fi]))) # or use "list(zip(r, Gr)))"


xpdfsuite_q = test_xpdfsuite_qsq[:, 0]
xpdfsuite_sq = test_xpdfsuite_qsq[:, 1]
xpdfsuite_q = test_xpdfsuite_qfq[:, 0]
xpdfsuite_fq = test_xpdfsuite_qfq[:, 1]
xpdfsuite_r = test_xpdfsuite_rgr[:, 0]
xpdfsuite_gr = test_xpdfsuite_rgr[:, 1]


exp_qiq = np.loadtxt(expqIq_data, skiprows=4)
exp_q = exp_qiq[:, 0]
# print('exp_q = ', exp_q)
# print('length of exp_q = ', len(exp_q))
exp_Iq = exp_qiq[:, 1]

plt.figure(1)
plt.plot(q, (scaled_expIq), label='exp_Iq')
plt.plot(q, (list_scaled_bkgIq), label='nn1 * bkg')
plt.plot(q, (Iq), label='bkg sbutracted exp_Iq')
plt.plot(list_q, (list_compton_scat), label='compton_scat')
plt.plot(list_q, (mean_sq_fi), label='mean_square_fi')
plt.plot(list_q, (sq_mean_fi), label='square_mean_fi')
plt.plot(list_q, 12*(list_compton_scat + mean_sq_fi), label='compton_scat + mean_sq_fi')
plt.xlabel('q (1/A)')
plt.ylabel('I(q)')
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(q, list_Sq, label='not normalized S(q)')
plt.plot(xpdfsuite_q, xpdfsuite_sq, label='xpdfsuite_S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.grid()
plt.legend()

plt.figure(3)
plt.plot(q, Sq, label='S(q)')
plt.plot(xpdfsuite_q, xpdfsuite_sq, label='xpdfsuite_S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.grid()
plt.legend()

plt.figure(4)
plt.plot(q, Fq, label='F(q)')
plt.plot(xpdfsuite_q, xpdfsuite_fq, label='xpdfsuite_F(q)')
plt.xlabel('q (1/A)')
plt.ylabel('F(q)')
plt.grid()
plt.legend()


plt.figure(5)
plt.plot(q, Iq/12, label='bkg sbutracted Iq')
plt.plot(q, mean_sq_fi, label='<f^2>')
plt.plot(q, sq_mean_fi, label='<f>^2')
plt.xlabel('q (1/A)')
plt.ylabel('<f^2> <f>^2')
plt.grid()
plt.legend()

#get G(r) data from integral and IFFT methods
r1, Gr1 = proc.calculate_expGr_integral(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep)
r2, Gr2 = proc.calculate_expGr_fft(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep, extrapolate_type="linear")

plt.figure(6)
plt.plot(r1, Gr1, label='integral-G(r)')
plt.plot(r2, Gr2, label='ifft-G(r)')
plt.plot(xpdfsuite_r, xpdfsuite_gr, label='xpdfsuite_G(r)' )
plt.xlabel('r (A)')
plt.ylabel('G(r)')
plt.grid()
plt.legend()

plt.figure(7)
plt.plot(xpdfsuite_r, xpdfsuite_gr, label='xpdfsuite_G(r)' )
plt.xlabel('r (A)')
plt.ylabel('G(r)')
plt.grid()
plt.legend()

plt.figure(8)
plt.plot(xpdfsuite_q, xpdfsuite_sq, label='xpdfsuite_S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.grid()
plt.legend()

plt.figure(9)
plt.plot(xpdfsuite_q, xpdfsuite_fq, label='xpdfsuite_F(q)')
plt.xlabel('q (1/A)')
plt.ylabel('F(q)')
plt.grid()
plt.legend()
plt.show()