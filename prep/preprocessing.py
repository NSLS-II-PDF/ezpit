import numpy as np
from scipy.spatial.distance import cdist
from scipy import interpolate


def create_atom_distance_matrix(atom_positions):
    distance_matrix = cdist(atom_positions, atom_positions)
    return distance_matrix


def __calc_fi(scat_values, q):
    """
    Support function
    """
    k_sq = (0.25 * q / np.pi) ** 2
    fi1 = scat_values[0] * np.exp(-scat_values[1] * k_sq)
    fi2 = scat_values[2] * np.exp(-scat_values[3] * k_sq)
    fi3 = scat_values[4] * np.exp(-scat_values[5] * k_sq)
    fi4 = scat_values[6] * np.exp(-scat_values[7] * k_sq)
    fic = scat_values[8]
    fi = fi1 + fi2 + fi3 + fi4 + fic
    return fi


def calculate_Iq(atom_indices, scattering_factors, atom_distance_matrix,
                 qmin=0.5, qmax=20, qstep=0.05):
    num_atom = len(atom_indices)
    num_fact = len(scattering_factors)
    diag_idx = np.diag_indices(num_atom)
    distance_matrix_non_zero = np.copy(atom_distance_matrix)
    distance_matrix_non_zero[diag_idx] = 1.0
    list_Iq = []
    q_range = np.arange(qmin, qmax, qstep)
    for q in q_range:
        fi_mat = np.zeros((num_atom, num_atom))
        list_fi = []
        for k in range(num_fact):
            list_fi.append(__calc_fi(scattering_factors[k], q))
        list_fi = np.asarray(list_fi)
        for i, idx in enumerate(atom_indices):
            fi = list_fi[idx]
            fi_mat[i, :] = fi
        if q == 0:
            sin_mat = np.ones(np.shape(atom_distance_matrix))
        else:
            sin_mat = np.sin(q * atom_distance_matrix) / (
                    q * distance_matrix_non_zero)
        sin_mat[diag_idx] = 1.0
        Iq = fi_mat * np.transpose(fi_mat) * sin_mat
        list_Iq.append(np.sum(Iq))
    return q_range, np.asarray(list_Iq)


def calculate_Sq(atom_indices, scattering_factors, atom_distance_matrix,
                 qmin=0.5, qmax=20, qstep=0.05, return_Iq=False):
    num_atom = len(atom_indices)
    num_fact = len(scattering_factors)
    print('num_fact = ', num_fact)
    diag_idx = np.diag_indices(num_atom)
    print('diag_idx = ', diag_idx)
    distance_matrix_non_zero = np.copy(atom_distance_matrix)
    distance_matrix_non_zero[diag_idx] = 1.0

    list_Iq = []
    sq_mean_fi = []
    mean_sq_fi = []
    q_range = np.arange(qmin, qmax, qstep)
    for q in q_range:
        fi_mat = np.zeros((num_atom, num_atom))
        list_fi = []
        for k in range(num_fact):
            list_fi.append(__calc_fi(scattering_factors[k], q))
        list_fi = np.asarray(list_fi)
        fi_sum, fi2_sum = 0.0, 0.0
        for i, idx in enumerate(atom_indices):
            fi = list_fi[idx]
            fi_mat[i, :] = fi
            fi_sum = fi_sum + fi
            fi2_sum = fi2_sum + fi ** 2
        if q == 0:
            sin_mat = np.ones(np.shape(atom_distance_matrix))
        else:
            sin_mat = np.sin(q * atom_distance_matrix) / (
                    q * distance_matrix_non_zero)
        sin_mat[diag_idx] = 1.0
        Iq = fi_mat * np.transpose(fi_mat) * sin_mat
        list_Iq.append(np.sum(Iq))
        sq_mean_fi.append((fi_sum / num_atom) ** 2)
        mean_sq_fi.append(fi2_sum / num_atom)
    list_Iq = np.asarray(list_Iq)
    mean_sq_fi = np.asarray(mean_sq_fi)  # <f^2>
    sq_mean_fi = np.asarray(sq_mean_fi)  # <f>^2
    list_Sq = (list_Iq - num_atom * mean_sq_fi) / (num_atom * sq_mean_fi) + 1
    list_Fq = q_range * (list_Sq - 1)  # added 8/26/2023

    if return_Iq:
        return q_range, list_Iq, list_Sq, list_Fq, mean_sq_fi, sq_mean_fi
    else:
        return q_range, list_Sq


def calculate_Gr_integral(q, Sq, rmin=0, rmax=100, rstep=0.02, qdamp=0.0):
    list_r = np.arange(rmin, rmax + rstep, rstep)
    Fq = (Sq - 1) * q
    qstep = q[1] - q[0]
    list_Gr = np.zeros_like(list_r)
    list_qFq = Fq * qstep
    for i, r in enumerate(list_r):
        list_Gr[i] = np.sum(list_qFq * np.sin(q * r)) * 2.0 / np.pi
    if qdamp != 0.0:
        list_Gr = np.exp(-0.5 * (list_r * qdamp) ** 2) * list_Gr
        print('list_gr with qdamp = ', list_Gr)
    return list_r, list_Gr


def calculate_Gr_fft(q, Sq, rmin=0, rmax=100, rstep=0.02, qdamp=0.0,
                     extrapolate_type="linear"):
    qstep = q[1] - q[0]
    print('q[0] = ', q[0])
    num_point = int(np.ceil(q[0] / qstep))
    print('num_point = ', num_point)
    Fq = (Sq - 1) * q
    pad_Fq = np.zeros(num_point)
    print('pad_Fq = ', pad_Fq)
    print('lenght of pad_Fq = ', len(pad_Fq))
    print('Fq before pad = ', Fq)
    print('lenght of Fq before pad = ', len(Fq))
    if q[0] > 0.0:
        f_inter = interpolate.interp1d(q, Fq, fill_value="extrapolate",
                                       kind=extrapolate_type)
        pad_Fq[:num_point] = f_inter(np.arange(num_point) * qstep)
    print('lenght of pad_Fq after interpolation = ', len(pad_Fq))
    Fq_vals = np.append(pad_Fq, Fq)
    print('Fq_vals = ', Fq_vals)
    print('lenght of Fq_vals = ', len(Fq_vals))
    r_list = np.arange(rmin, rmax + rstep, rstep)
    num_point = len(Fq_vals)
    print('num_point of len(Fq_vals) = ', num_point)
    total_point = int(2 * np.pi / (rstep * qstep))
    print('total_point = ', total_point)
    Fq_pad = np.pad(Fq_vals, (0, total_point - num_point), mode="constant")
    print('length of Fq_pad = ', Fq_pad)
    norm = total_point * qstep * 2 / np.pi
    gr_fine = norm * np.imag(np.fft.ifft(Fq_pad))
    rfine = np.arange(total_point) * rstep
    if qdamp != 0.0:
        gr_fine = gr_fine * np.exp(-0.5 * (rfine * qdamp) ** 2)
    gr = np.interp(r_list, rfine, gr_fine)
    return r_list, gr


def cost_function(calgr_norm, expgr_norm):
    return np.sqrt(sum((calgr_norm - expgr_norm) ** 2) / sum((expgr_norm) ** 2))
