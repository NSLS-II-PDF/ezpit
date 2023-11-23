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
    #print('num_fact = ', num_fact)
    diag_idx = np.diag_indices(num_atom)
    #print('diag_idx = ', diag_idx)
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
        #print('list_gr with qdamp = ', list_Gr)
    return list_r, list_Gr


def calculate_Gr_fft(q, Sq, rmin=0, rmax=100, rstep=0.02, qdamp=0.0,
                     extrapolate_type="linear"):
    qstep = q[1] - q[0]
    #print('q[0] = ', q[0])
    num_point = int(np.ceil(q[0] / qstep))
    #print('num_point = ', num_point)
    Fq = (Sq - 1) * q
    pad_Fq = np.zeros(num_point)
    #print('pad_Fq = ', pad_Fq)
    #print('lenght of pad_Fq = ', len(pad_Fq))
    #print('Fq before pad = ', Fq)
    #print('lenght of Fq before pad = ', len(Fq))
    if q[0] > 0.0:
        f_inter = interpolate.interp1d(q, Fq, fill_value="extrapolate",
                                       kind=extrapolate_type)
        pad_Fq[:num_point] = f_inter(np.arange(num_point) * qstep)
    #print('lenght of pad_Fq after interpolation = ', len(pad_Fq))
    Fq_vals = np.append(pad_Fq, Fq)
    #print('Fq_vals = ', Fq_vals)
    #print('lenght of Fq_vals = ', len(Fq_vals))
    r_list = np.arange(rmin, rmax + rstep, rstep)
    num_point = len(Fq_vals)
    #print('num_point of len(Fq_vals) = ', num_point)
    total_point = int(2 * np.pi / (rstep * qstep))
    #print('total_point = ', total_point)
    Fq_pad = np.pad(Fq_vals, (0, total_point - num_point), mode="constant")
    #print('length of Fq_pad = ', Fq_pad)
    norm = total_point * qstep * 2 / np.pi
    gr_fine = norm * np.imag(np.fft.ifft(Fq_pad))
    rfine = np.arange(total_point) * rstep
    if qdamp != 0.0:
        gr_fine = gr_fine * np.exp(-0.5 * (rfine * qdamp) ** 2)
    gr = np.interp(r_list, rfine, gr_fine)
    return r_list, gr

#calculate S(q), F(q), G(r) from experimental I(q)
def calculate_expSq(atom_indices, scattering_factors, expqIq_data, bkgqIq_data, qmin=0, qmax=25, qstep=0.1,
                    background_scale=1.1, poly_order=11, return_Iq=False):
    # load experimental Iq data
    exp_qiq = np.loadtxt(expqIq_data)
    exp_q = exp_qiq[:, 0]
    exp_Iq = exp_qiq[:, 1]

    #load background data   #####think about adding more background data
    bkg_qiq = np.loadtxt(bkgqIq_data)
    bkg_Iq = bkg_qiq[:, 1]

    num_atom = len(atom_indices)
    num_fact = len(scattering_factors)

    sq_mean_fi = []
    mean_sq_fi = []
    q_range = np.arange(qmin, qmax, qstep)

    #scaled expdata and bkgdata before subtract each other################11/02/2023
    scaled_expIq = np.interp(q_range, exp_q, exp_Iq)
    scaled_bkgIq = np.interp(q_range, exp_q, bkg_Iq)
    list_scaled_bkgIq = background_scale * scaled_bkgIq
    #bkg subtracted data
    list_Iq = scaled_expIq - list_scaled_bkgIq

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
            fi2_sum = fi2_sum + (fi ** 2)
        sq_mean_fi.append((fi_sum / num_atom) ** 2)
        mean_sq_fi.append(fi2_sum / num_atom)
    #list_Iq = np.asarray(list_Iq)
    mean_sq_fi = np.asarray(mean_sq_fi)  # <f^2>
    sq_mean_fi = np.asarray(sq_mean_fi)  # <f>^2
    list_Sq = (list_Iq - num_atom * mean_sq_fi) / (num_atom * sq_mean_fi) + 1
    #list_Sq = (list_Iq - mean_sq_fi) / (sq_mean_fi) + 1

    #calculate_norm_expSq(q_range, list_Sq, poly_order = 11):
    # generate two numpy array [[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    #                         [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]
    # poly_fit = np.zeros((2, poly_order), dtype=float) ---> commented out on Aug11, 2023
    # np.polyfit()
    qmax = np.max(q_range)
    qmax_for_sq = qmax
    #print('qmax = ', qmax)
    #print('qmax_for_sq = ', qmax_for_sq)
    rpoly = np.pi * poly_order / qmax_for_sq
    #print('rpoly = ', rpoly)

    qscaled_by_qmax_for_sq = q_range / qmax_for_sq  # qscaled_by_qmax_for_sq is numpy array
    matrix_vandermonde = np.vander(qscaled_by_qmax_for_sq, poly_order + 1)
    # remove 1 (the last element =a0) in vandermonde matrix:
    # 4.04E-17 1.25E-15 3.86E-14 1.19E-12 3.70E-11 1.14E-09 3.53E-08 1.09E-06 3.38E-05 1.05E-03 3.23E-02
    remove_a0_matrix_vandermonde = matrix_vandermonde[:, :-1]
    # remove_a0_matrix_vandermonde = matrix_vandermonde
    # fit a polynomial, minimizing least square difference using  lstsq_fq = q * (sq_from_sq_dataframe - 1.0)
    # ***********lstsq_fq is fq obtained from unnormalized S(q)************
    # ***********bring unnormalized fq osciallation is going to zero--> thus make 1 difference between S(q) and F(q)
    lstsq_fq = q_range * (list_Sq - 1.0)  # experimenatl F(q) = Fm(q)
    # p: ndarray, least-squares solution, residuals: ndarray, Sums of squared residuals, rank: int, Rank of matrix a
    # s: ndarray, Singular values of a

    # https://jyyuan.wordpress.com/2014/01/02/polynomial-interpolation-using-vandermonde-matrix-and-least-squares/
    # vandermonde matrix works with "np.linalg.lstsq" to minimize least-sqare difference. August/15/2023.
    # So aside from the downsides of interpolation and whatever, what should we be careful of here?
    # Overfitting the polynomial can make for some very poor solutions that don’t really make any sense in the context of the problem at hand,
    # so in general, doing least squares with a tall Vandermonde matrix for this interpolation problem will get better results than a square Vandermonde or an underdetermined problem.
    # The saying that simpler is better has never rung so true.
    # A more sensitive issue specific to computer calculation is that the Vandermonde matrix built can have some poor conditioning,
    # with larger powers ballooning or shrinking and causing instability of solutions. In those cases,
    # an alternative and numerical stable interpolation method using Lagrange interpolation is preferred.
    p, residuals, rank, s = np.linalg.lstsq(remove_a0_matrix_vandermonde, lstsq_fq, rcond=None)

    # vanderMat = np.vander([5],5), print('vanderMat = ', vanderMat) --> vanderMat =  [[625 125  25   5   1]]
    # generate one-dimensional matrix --> thus can compute p (one-dimension) / vandermonde_before_2nd_vander (ond-dimension)
    vandermonde_2nd = (np.vander([qmax_for_sq], poly_order + 1)[0, :-1])  # -1 plays a role to remove unit (1) in vandermonde matrix
    new_p = p / vandermonde_2nd

    # poly1d just add x^n in the coefficient obtained from the vandermonde matrix
    # p_after_vander =  [-4.37409054e-13, 6.60397573e-11, -4.37727463e-09]
    # --> final_p = -4.374e-13 x^2 + 6.604e-11 x - 4.377e-09
    final_p = np.poly1d(new_p)
    # get polynomial data as function of q
    polynomial_for_sq = final_p(q_range)
    norm_list_Sq = list_Sq - polynomial_for_sq
    list_Fq = q_range * (norm_list_Sq - 1)

    if return_Iq:
        return q_range, list_Iq, scaled_expIq, list_scaled_bkgIq, list_Sq, norm_list_Sq, list_Fq, mean_sq_fi, sq_mean_fi
    else:
        return q_range, list_Iq, scaled_expIq, list_scaled_bkgIq, list_Sq, norm_list_Sq, list_Fq, mean_sq_fi, sq_mean_fi


#calculate_expGr_integral has no function of qdamp
def calculate_expGr_integral(q, list_Sq, rmin=0, rmax=100, rstep=0.02):
    list_r = np.arange(rmin, rmax + rstep, rstep)
    Fq = (list_Sq - 1) * q
    qstep = q[1] - q[0]
    list_Gr = np.zeros_like(list_r)
    list_qFq = Fq * qstep
    for i, r in enumerate(list_r):
        list_Gr[i] = np.sum(list_qFq * np.sin(q * r)) * 2.0 / np.pi
    return list_r, list_Gr


def calculate_expGr_fft(q, Sq, rmin=0, rmax =100, rstep=0.01,
                        extrapolate_type="linear"): #"extrapolation" and "zeropad"
    qstep = q[1] - q[0]
    #print('q[0] = ', q[0])
    num_point = int(np.ceil(q[0] / qstep))
    #print('num_point = ', num_point)
    Fq = (Sq - 1) * q
    pad_Fq = np.zeros(num_point)
    #print('pad_Fq = ', pad_Fq)
    #print ('lenght of pad_Fq = ', len(pad_Fq))
    #print('Fq before pad = ', Fq)
    #print('lenght of Fq before pad = ', len(Fq))

    if q[0] > 0.0:
        f_inter = interpolate.interp1d(q, Fq, fill_value="extrapolate",
                                       kind=extrapolate_type)
        pad_Fq[:num_point] = f_inter(np.arange(num_point) * qstep)
    #print('lenght of pad_Fq after interpolation = ', len(pad_Fq))
    Fq_vals = np.append(pad_Fq, Fq)
    #print('Fq_vals = ', Fq_vals)
    #print('lenght of Fq_vals = ', len(Fq_vals))
    r_list = np.arange(rmin, rmax + rstep, rstep)
    num_point = len(Fq_vals)
    #print('num_point of len(Fq_vals) = ', num_point)
    total_point = int(2 * np.pi / (rstep * qstep))
    #print('total_point = ', total_point)
    Fq_pad = np.pad(Fq_vals, (0, total_point - num_point), mode="constant")
    #print('length of Fq_pad = ', Fq_pad)
    norm = total_point * qstep * 2 / np.pi
    gr_fine = norm * np.imag(np.fft.ifft(Fq_pad))
    rfine = np.arange(total_point) * rstep
    gr = np.interp(r_list, rfine, gr_fine)
    return r_list, gr

def cost_function(calgr_norm, expgr_norm):
    return np.sqrt(sum((calgr_norm - expgr_norm) ** 2) / sum((expgr_norm) ** 2))
