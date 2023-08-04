import numpy as np
import spectral_modelling.utils.constants as _c
from spectral_modelling.model_and_invert.fas_log import ln_fas
import spectral_modelling.utils.utils as utils


def handle_func_deriv(p_in, *fcn_args):
    """
    Calculate the jacobian vector of the handle function (with possibility to 
    use custom weights and scaled for total number of points)

    NOTE All quantities are in array format (no matrices), 
         the output is an array
    """

    stat_p = fcn_args[0]
    index_v = fcn_args[1]
    p_fixed = fcn_args[2]
    index_f = fcn_args[3]
    N_gamma = fcn_args[4]
    N_delta = fcn_args[5]
    utils.checkext_statparams(stat_p)

    M = stat_p.M
    F = stat_p.F
    R = stat_p.R
    N_ev = stat_p.N_ev
    N_sta = stat_p.N_sta
    N_freq = stat_p.N_freq

    if not (isinstance(p_in, np.ndarray)):
        p_in = p_in.pars

    # rebuild total parameter array and extract some parameter blocks
    p_input = utils.insert_elem(p_in, index_f, p_fixed)

    beta = p_input[N_ev:N_ev + N_ev]
    beta = np.tile(np.array([beta]).T, (1, N_sta * N_freq))
    beta = beta.reshape(-1)
    gamma = p_input[(N_ev + N_ev):(N_ev + N_ev + N_gamma)]
    delta = p_input[(N_ev + N_ev + N_gamma):(N_ev + N_ev + N_gamma + N_delta)]
    kappa = p_input[
            (N_ev + N_ev + N_gamma + N_delta + N_sta):(N_ev + N_ev + N_gamma
                                                       + N_delta + 2 * N_sta)]
    kappa = np.tile(np.array([kappa]).T, N_freq).reshape(-1)
    kappa = np.tile(kappa, N_ev)

    wresd = np.zeros_like(p_input)
    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    diff = M * ln_fas(p_input, fcn_args_ee, diff=True)

    weights = stat_p.weights
    if not isinstance(weights, np.ndarray):
        weights = np.ones_like(diff)
    if len(weights) != len(diff):
        raise utils.MinimizeException("Weights array must be of the same lenght"
                                      " as data")

    # calculate jacobian blocks for each corresponding parameter block
    # aplha
    diff_alpha = diff*weights
    for l in range(N_ev):
        wresd[l] = 2. * (diff_alpha[N_sta*N_freq*l:N_sta*N_freq*(l + 1)].sum())

    # beta
    diff_beta = diff * (F ** 2 / (beta*(beta ** 2 + F ** 2))) * weights
    for l in range(N_ev):
        wresd[N_ev + l] = 4. * (diff_beta[N_sta*N_freq*l:
                                          N_sta*N_freq*(l + 1)].sum())

    # gamma
    diff_gamma = diff*weights
    for l in range(N_gamma):
        diff_gamma_l = diff_gamma * utils.ln_piecew_func(R, F, gamma, der=l)
        wresd[N_ev + N_ev + l] = 2. * diff_gamma_l.sum()

    # delta
    diff_delta0 = diff * (F ** (1. - delta[1])) * R * weights
    wresd[N_ev + N_ev + N_gamma] = (2.*_c.PI / (_c.V_S*(delta[0] ** 2))) \
                                   * diff_delta0.sum()

    diff_delta1 = diff * (F ** (1. - delta[1])) * np.log(F) \
                  * (R / (_c.V_S*delta[0]) + kappa) * weights
    wresd[N_ev + N_ev + N_gamma + 1] = 2. * _c.PI * diff_delta1.sum()

    # site_fi
    diff_site_fi = diff * weights
    for l in range(N_sta):
        dsf = 0
        for i in range(N_ev):
            for k in range(N_freq):
                dsf += (diff_site_fi[N_sta*N_freq*i + N_freq*l + k])
        wresd[N_ev + N_ev + N_gamma + N_delta + l] = 2.*dsf

    # site_k
    diff_site_k = diff * (F ** (1. - delta[1])) * weights
    for l in range(N_sta):
        dsk = 0
        for i in range(N_ev):
            for k in range(N_freq):
                dsk += diff_site_k[N_sta*N_freq*i + N_freq*l + k]
        wresd[N_ev + N_ev + N_gamma + N_delta + N_sta + l] = -(_c.PI*2.) * dsk

    # epsilons
    diff_eps = diff*weights
    for l in range(N_ev):
        wresd[N_ev + N_ev + N_gamma + N_delta + 2*N_sta + l] = \
            2. * (diff_eps[N_sta*N_freq*l:N_sta*N_freq*(l + 1)].sum())

    wresd[N_ev + N_ev + N_gamma + N_delta + 2*N_sta + N_ev] = 2.*diff_eps.sum()

    for l in range(N_sta):
        esi = 0
        for i in range(N_ev):
            for k in range(N_freq):
                esi += diff_eps[N_sta*N_freq*i + N_freq*l + k]
        wresd[
            N_ev + N_ev + N_gamma + N_delta + 2*N_sta + N_ev + 1 + l] = 2.*esi

    # scale and keep components related to free parameters only
    wresd = wresd / (N_ev * N_sta * N_freq)
    wresd = wresd[index_v]
    return wresd
