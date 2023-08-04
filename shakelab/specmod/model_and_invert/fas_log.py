# =============================================================================
#
# Place for Copyright and authors information
#
# =============================================================================
"""
Implementation of the forward modelling for Fourier 
amplitude spectra for seismic traces, based on the formulas by
Bora et al. (2017), Franceschina et al. (2006) and Edwards et al. (2008)

cf. Equation (4-12) of PhD thesis (Cataldi)

NB: frequency-dependent site contributions a(f) are not modeled and
left to be calculated from the residuals
"""

import numpy as np
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils

def scale_to_motion(motion, F):
    if motion.lower() == 'displacement':
        return 0
    elif motion.lower() == 'velocity':
        return np.log(2*_c.PI*F)
    elif motion.lower() == 'acceleration':
        return 2*np.log(2*_c.PI*F)
    else:
        return None


# =============================================================================

def ln_alpha(alpha_in, N_sta, N_freq):
    """
    Compute the long-period plateau log-value of the source spectrum at
    the source (log(M0))

    :param np.array(N_ev) alpha_in:
        Array, input log(M0) values for each event

    :param int N_sta, N_freq:
        The number of stations and frequency points in the dataset

    :return np.array(N_ev*N_sta*N_freq) alpha:
        Array, log of the source momentum M0
    """

    alpha = np.tile(np.array([alpha_in]).T, (1, N_sta*N_freq))

    alpha = alpha.reshape(-1)

    return alpha


# =============================================================================

def ln_source_spectrum(beta_in, F, N_sta, N_freq):
    """
    Compute the source spectrum shape using the Brune omega-square
    model

    :param np.array(N_ev) beta_in:
        Array, input corner frequency (f_c) values for each event

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param int N_sta, N_freq:
        The number of stations and frequency points in the dataset

    :return np.array(N_ev*N_sta*N_freq) ss:
        Array, log of Brune source spectrum shape
    """

    beta = np.tile(np.array([beta_in]).T, (1, N_sta*N_freq))

    beta = beta.reshape(-1)

    ss = -np.log(1 + (F/beta)**2)

    return ss


# =============================================================================

def ln_app_geom_spreading(gamma_in, F, R):
    """
    Compute the apparent geometrical spreading piecewise function
    (from Franceschina et al., 2006)

    :param np.array(N_gamma) gamma_in:
        Array, input exponent values

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param np.array(N_ev*N_sta*N_freq) R:
        Array, used hypocentral distances

    :return np.array(N_ev*N_sta*N_freq) ags:
        Array, log of apparent geometrical spreading
    """

    ags = utils.ln_piecew_func(R, F, gamma_in, der=None)

    return ags


# =============================================================================

def ln_q_attenuation(delta_in, F, R):
    """
    Compute Q dependent attenuation contribution to Fourier 
    amplitude spectrum

    :param int delta_in:
        The input [Q0, Q_alpha] value

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param np.array(N_ev*N_sta*N_freq) R:
        Array, used hypocentral distances

    :return np.array(N_ev*N_sta*N_freq) q_a:
        Array, log of Q-dependent attenuation function
    """

    Q_0 = delta_in[0]
    Q_eta = delta_in[1]
    q_att = -(_c.QCOST/Q_0)*(F**(1.-Q_eta))*R

    return q_att


# =============================================================================

def ln_source(alpha_in, beta_in, F, N_sta, N_freq):
    """
    Compute the collective signal moment

    :param np.array(N_ev) alpha_in:
        Array, input log(M0) values for each event

    :param np.array(N_ev) beta_in:
        Array, input corner frequency (f_c) values for each event

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param int N_sta, N_freq:
        The number of stations and frequency points in the dataset

    :return np.array(N_ev*N_sta*N_freq) soc:
        Array, collective source contribution to FAS
    """

    soc = (_c.LN_MCOST + ln_alpha(alpha_in, N_sta, N_freq) 
           + ln_source_spectrum(beta_in, F, N_sta, N_freq))

    return soc

# =============================================================================


def ln_path(gamma_in, delta_in, F, R):
    """
    Compute the collective signal moment

    :param np.array(N_gamma) gamma_in:
        Array, input exponent values

    :param int delta_in:
        The input Q0 val

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param np.array(N_ev*N_sta*N_freq) R:
        Array, used hypocentral distances

    :return np.array(N_ev*N_sta*N_freq) pac:
        Array, collective path contribution to FAS
    """

    pac = ln_app_geom_spreading(gamma_in, F, R) \
          + ln_q_attenuation(delta_in, F, R)

    return pac


# =============================================================================

def ln_site(site_fi_in, site_k_in, delta_in, F, N_ev, N_freq):
    """
    Compute the collective signal moment

    :param np.array(N_sta) site_fi_in:
        Array, input freq.-indep. amplification values for each site

    :param np.array(N_sta) site_k_in:
        Array, input k values for each site

    :param int delta_in:
        The input [Q0, Q_alpha] value

    :param np.array(N_ev*N_sta*N_freq) F:
        Array, used frequency points

    :param int N_ev, N_sta, N_freq:
        The number of events, stations and frequency points in the
        dataset

    :return np.array(N_ev*N_sta*N_freq) sic:
        The collective site contribution to FAS
    """

    Q_eta = delta_in[1]

    SFI = np.tile(np.array([site_fi_in]).T, N_freq).reshape(-1)
    SFI = np.tile(SFI, N_ev)

    SK = np.tile(np.array([site_k_in]).T, N_freq).reshape(-1)
    SK = np.tile(SK, N_ev)
    kappa = - _c.PI*(F**(1.-Q_eta))*SK

    sic = SFI + kappa

    return sic


# =============================================================================

def epsilons(eps_source_input, eps_path_input, eps_site_input, N_ev, N_sta,
             N_freq):
    """
    Add the uncertainty terms

    :param np.array(N_ev) eps_source_input:
        Array, uncertainty on source contribution

    :param np.array(1) eps_path_input:
        Single element array, uncertainty on path contribution

    :param np.array(N_sta) eps_site_input:
        Array, uncertainty on site contribution

    :param int N_ev, N_sta, N_freq:
        The number of events, stations and frequency points in the
        dataset

    :return np.array(N_ev*N_sta*N_freq) epsilons:
        The collective uncertainty to FAS
    """

    ESO = np.tile(np.array([eps_source_input]).T, N_sta*N_freq).reshape(-1)

    EPA = np.tile(np.array([eps_path_input]), N_ev*N_sta*N_freq).reshape(-1)

    ESI = np.tile(np.array([eps_site_input]).T, N_freq).reshape(-1)
    ESI = np.tile(ESI, N_ev)

    errors = ESO + EPA + ESI

    return errors


# =============================================================================

MOTION = cfg.MODELS["MOTION"]


def ln_fas(p_in, fcn_args, diff=False):
    """
    Calculate the forward model of the ln of Fourier amplitude spectrum 
    for each given event and site.
    First computes the theoretical (synthetic) spectra, regardless of the
    presence of data already stored in stat_p.data, then subtracts the 
    observed spectra if parameter 'diff' is True.
    If the scaling factor stat_p.scale is not null, the result is
    scaled.
    This version supports parameter fixing; information on which
    parameters are fixed and which are free to vary is included in
    fcn_args components.

    NOTE input data is transformed in log form for calculations

    Used parameters are:
          alpha <--> ln(M_0) [g * cm^2  / s^2]
          beta  <--> f_c [Hz]
          gamma <--> gamma indexes of geom. spread. function
                     [adimensional]
          delta <--> Q_0 [adimensional] (if Q_alpha is not used)
          site_fi <--> ln(A) [adimensional]
          site_k  <--> k attenuation factor [s]

    :return np.array(N_ev*N_sta*N_freq) fas:
        The fourier amplitude spectrum values for given F and R
    """

    stat_p, index_v, p_fixed, index_f, N_gamma, N_delta = fcn_args
    utils.checkext_statparams(stat_p)

    F = stat_p.F
    R = stat_p.R
    N_ev = stat_p.N_ev
    N_sta = stat_p.N_sta
    N_freq = stat_p.N_freq
    lndata = np.log(stat_p.data)
    scale = stat_p.scale

    if not (isinstance(p_in, np.ndarray)):
        p_in = p_in.pars

    p_input = utils.insert_elem(p_in, index_f, p_fixed)
    p_input = utils.array_to_Param(p_input, stat_p, N_gamma, N_delta)

    fas = (scale_to_motion(MOTION, F) 
           + ln_source(p_input.alpha, p_input.beta, F, N_sta, N_freq)
           + ln_path(p_input.gamma, p_input.delta, F, R)
           + ln_site(p_input.site_fi, p_input.site_k, p_input.delta, F, N_ev,
                     N_freq)
           + epsilons(p_input.eps_source, p_input.eps_path, p_input.eps_site,
                      N_ev, N_sta, N_freq))

    if scale is not None:
        fas = fas / scale
        lndata = lndata / scale
    
    if diff is False:
        return fas
    else:
        return fas - lndata


def handle_func(p_in, *fcn_args):
    """
    Calculate the handle (=cost) function as: 
        sum( ((Z_calc-data)**2)*weights )/(N_ev*N_sta*N_freq)

    NOTE All quantities are in array format (no matrices), 
         the output is a scalar
    """

    stat_p = fcn_args[0]
    index_v = fcn_args[1]
    p_fixed = fcn_args[2]
    index_f = fcn_args[3]
    N_gamma = fcn_args[4]
    N_delta = fcn_args[5]
    utils.checkext_statparams(stat_p)

    fcn_args = (stat_p, index_v, p_fixed, index_f, N_gamma, N_delta)

    diff = stat_p.M * ln_fas(p_in, fcn_args, diff=True)
    weights = stat_p.weights
    # if no weights are used, set the weight array to unit array
    if not isinstance(weights, np.ndarray):
        weights = np.ones_like(diff)

    if len(weights) != len(diff):
        raise utils.MinimizeException("Weights array must be of the same lenght"
                                      " as data")

    cost = (np.square(diff)*weights).sum()
    cost = cost/(stat_p.N_ev*stat_p.N_sta*stat_p.N_freq)

    return cost
