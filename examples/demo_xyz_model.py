import matplotlib.pyplot as plt
import losa.loadersaver as losa
import proc.processing as proc
import timeit

input_base = "C:/Users/gkwon/Pycharmprojects/ezpit/data/"

input_exp_file = input_base + '/D_S841-00003.gr'
atom_xyz_file = input_base + "/5IrC_r5a-1Ir.xyz" #5IrC_r5a-1Ir.xyz"  #Ni(OH)2-109391-ICSD-10x10x1.xyz"
aff_element_file = input_base + '/aff_elementonly.txt'
aff_parm_file = input_base + 'aff_parmonly.txt'

atom_names, atom_positions = losa.load_atom_name_positions(atom_xyz_file)
#print(atom_names)

#import sys
#sys.exit(0)
# Database
database_atom_names = losa.load_atom_names(aff_element_file)
database_scat_factors = losa.load_scattering_factors(aff_parm_file)

# print(atom_names)
# print(atom_positions)
# print(len(database_atom_names))
# print(database_scat_factors.shape)

qmin = 0.0
qmax = 24
qstep = 0.01
rmin = 0
rmax = 12
rstep = 0.001
qdamp = 0.0
#method = 'ifft'

t5 = timeit.default_timer()
atom_distance_matrix = proc.create_atom_distance_matrix(atom_positions)
atom_unique_names, atom_counts, atom_indices = losa.group_atoms(atom_names)

scattering_factors = losa.get_scattering_factors(atom_unique_names, database_atom_names,
                                                database_scat_factors)

q, Iq, Sq, Fq, mean_sq_fi, sq_mean_fi = proc.calculate_Sq(atom_indices, scattering_factors, atom_distance_matrix,
                                    qmin=qmin, qmax=qmax, qstep=qstep, return_Iq=True)

# print(qIq.shape)
plt.figure(0)
plt.plot(q, Iq, label='I(q)')
#plt.plot(q, sq_mean_fi, label='sq_mean_fi')
#plt.plot(q, mean_sq_fi, label='mean_sq_fi')
plt.xlabel('q (1/A)')
plt.ylabel('I(q)')
plt.yscale("log")
plt.grid()
plt.legend()

plt.figure(1)
plt.plot(q, Sq, label='S(q)')
plt.xlabel('q (1/A)')
plt.ylabel('S(q)')
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(q, Fq, label='F(q)')
plt.xlabel('q (1/A)')
plt.ylabel('F(q)')
plt.grid()
plt.legend()

t0 = timeit.default_timer()
r, Gr = proc.calculate_Gr_integral(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep, qdamp=qdamp)
#np.savetxt(input_base + 'Gr_qdamp0p0.dat', np.column_stack(([r, Gr])))   # or use "list(zip(r, Gr)))"

t1 = timeit.default_timer()
print("Time cost real space!!! ", t1-t0)


plt.figure(3)
plt.plot(r, Gr, label='integral')
plt.xlabel('r (A)')
plt.ylabel('G(r)')
plt.legend()

t0 = timeit.default_timer()
r, Gr2 = proc.calculate_Gr_fft(q, Sq, rmin=rmin, rmax=rmax, rstep=rstep, qdamp=qdamp,
                              extrapolate_type="linear")
t1 = timeit.default_timer()
print("Time cost Fourier space!!! ", t1-t0)
t6 = timeit.default_timer()

print("Total time cost including real and Fourier space!!! ", t6-t5)

plt.plot(r, Gr2, label='ifft')
plt.xlabel('r (A)')
plt.ylabel('G(r)')
plt.grid()
plt.legend()

plt.figure(4)
plt.plot(r, Gr-Gr2, label='difference')
plt.xlabel('r (A)')
plt.ylabel('G(r)-G(r)2')
plt.grid()
plt.legend()
plt.show()
