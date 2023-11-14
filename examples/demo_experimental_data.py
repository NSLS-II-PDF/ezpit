import matplotlib.pyplot as plt
import losa.loadersaver as losa
import proc.processing as proc
import timeit
import numpy as np

input_base = "C:/Users/gkwon/PycharmProjects/ezpit/data/"
aff_element_file = input_base + 'aff_elementonly.txt'
aff_parm_file = input_base + 'aff_parmonly.txt'

# import sys
# sys.exit(0)
composition = {'Co':13, 'O':44, 'P': 1}
atom_names = losa.convert_atom_names(composition)
#load experimental Iq (expqIq_data and background I(q) (bkgqIq_data)
expqIq_data = 'C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5.chi'
bkgqIq_data = 'C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_0p7cap_Nsum6.chi'
#Parameters
qmin = 0.6
qmax = 24
qstep = 0.01
background_scale = 0.99
qdamp = 0
poly_order = 11
rmin = 0
rmax = 20
rstep = 0.01

method = 'ifft'

# Database
database_atom_names = losa.load_atom_names(aff_element_file)
database_scat_factors = losa.load_scattering_factors(aff_parm_file)

atom_unique_names, atom_counts, atom_indices = losa.group_atoms(atom_names)

scattering_factors = losa.get_scattering_factors(atom_unique_names, database_atom_names,
                                                database_scat_factors)
q, Iq, scaled_expIq, list_scaled_bkgIq, list_Sq, Sq, Fq, mean_sq_fi, sq_mean_fi = proc.calculate_expSq(atom_unique_names, atom_counts, atom_indices,
                                                             scattering_factors, expqIq_data, bkgqIq_data, qmin=qmin, qmax=qmax, qstep=qstep,
                                                             background_scale=background_scale, poly_order=11, return_Iq=False)
#load test results in xPDFsuite with the same data plotted here.
test_xpdfsuite_qsq = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-Nohdeader.sq")
test_xpdfsuite_qfq = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-NOheader.fq")
test_xpdfsuite_rgr = np.loadtxt("C:/Users/gkwon/Pycharmprojects/ezpit/data/sum_A_CoPiITO_110320-1_Nsum5-Noheader.gr")

xpdfsuite_q = test_xpdfsuite_qsq[:, 0]
xpdfsuite_sq = test_xpdfsuite_qsq[:, 1]
xpdfsuite_q = test_xpdfsuite_qfq[:, 0]
xpdfsuite_fq = test_xpdfsuite_qfq[:, 1]
xpdfsuite_r = test_xpdfsuite_rgr[:, 0]
xpdfsuite_gr = test_xpdfsuite_rgr[:, 1]

exp_qiq = np.loadtxt(expqIq_data)
exp_q = exp_qiq[:, 0]
# print('exp_q = ', exp_q)
# print('length of exp_q = ', len(exp_q))
exp_Iq = exp_qiq[:, 1]

plt.figure(1)
plt.plot(q, np.log(scaled_expIq), label='exp_Iq')
plt.plot(q, np.log(list_scaled_bkgIq), label='nn1 * bkg')
plt.plot(q, np.log(Iq), label='bkg sbutracted exp_Iq')
plt.xlabel('q (1/A)')
plt.ylabel('I(q)')
plt.legend()

plt.figure(2)
plt.plot(q, list_Sq, label='not normalized S(q)')
#plt.plot(xpdfsuite_q, xpdfsuite_sq, label='xpdfsuite_S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.legend()

plt.figure(3)
plt.plot(q, Sq, label='S(q)')
#plt.plot(xpdfsuite_q, xpdfsuite_sq, label='xpdfsuite_S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.legend()

plt.figure(4)
plt.plot(q, Fq, label='F(q)')
#plt.plot(xpdfsuite_q, xpdfsuite_fq, label='xpdfsuite_F(q)')
plt.xlabel('q (1/A)')
plt.ylabel('F(q)')
plt.legend()

plt.figure(5)
plt.plot(q, Iq/12, label='bkg sbutracted Iq')
plt.plot(q, mean_sq_fi, label='<f^2>')
plt.plot(q, sq_mean_fi, label='<f>^2')
plt.xlabel('q (1/A)')
plt.ylabel('<f^2> <f>^2')
plt.legend()

#get G(r) data from integral and IFFT methods
r1, Gr1 = proc.calculate_expGr_integral(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep)
r2, Gr2 = proc.calculate_expGr_fft(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep, extrapolate_type="linear")

plt.figure(6)
plt.plot(r1, Gr1, label='integral-G(r)')
plt.plot(r2, Gr2, label='ifft-G(r)')
#plt.plot(xpdfsuite_r, xpdfsuite_gr, label='xpdfsuite_G(r)' )
plt.xlabel('r (A)')
plt.ylabel('G(r)')
plt.legend()
plt.show()

