import numpy as np
import os
import pickle
import timeit
from copy import deepcopy
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils
import spectral_modelling.utils.myClass as myC
from spectral_modelling.model_and_invert.fas_log import handle_func
from spectral_modelling.model_and_invert.jacobian_module import handle_func_deriv
from spectral_modelling.model_and_invert.wrap_function import minimizer_wrapper

def bounds_NE_Italy(pars, alpha_ref, m0_ref):
    """
    Define initial parameter bounds for seismic setting of
    North East Italy
    """

    lb = []
    ub = []

    # alpha
    for i in range(len(pars.alpha)): 
        lb.append(alpha_ref[i] - 1.727)
        ub.append(alpha_ref[i] + 1.727)
    # beta
    for i in range(len(pars.beta)): 
        lb.append(utils.fcmin_NE_IT(m0_ref[i]))
        ub.append(utils.fcmax_NE_IT(m0_ref[i]))
    # gamma
    for i in range(len(pars.gamma)): 
        lb.append(-2.)
        ub.append(-0.3)
    # delta
    lb.append(20)
    ub.append(500)
    lb.append(0.)
    ub.append(1.)
    # site_fi
    for i in range(len(pars.site_fi)): 
        lb.append(-1.61)
        ub.append(1.61)
    # site_k
    for i in range(len(pars.site_k)): 
        lb.append(0.01)
        ub.append(0.15)
    # epsilons
    for i in range(len(pars.eps_source)+len(pars.eps_path)+len(pars.eps_site)): 
        lb.append(-np.inf)
        ub.append(np.inf)

    bounds_NE_IT = (lb, ub)
    return bounds_NE_IT


def fix_params(pars, model_name, use_uncert, ref_ampl, stas_dict):
    """
    Fix parameters in Params object according to chosen parameter 
    hypothesis (=model)
    """

    # initialize the 'vary' istance if not present
    if len(pars.vary) == 0:
        pars.vary_all(True)

    N_i = pars.alpha.size
    N_j = pars.site_fi.size
    N_gamma = pars.gamma.size
    N_delta = pars.delta.size

    # fixing Q0
    if model_name == 'malagnini_Q0':
        pars.set_vary((2*N_i+N_gamma+1,), False)

    # fixing all Q factors
    if model_name == 'gamma-deltas':
        for i in range(N_delta):
            pars.set_vary((2*N_i+N_gamma+i,), False)

    # case with fixed freq.-indep. amplification factors
    if model_name.startswith('fixA'):
        for j in range(N_j):
            pars.set_vary((2*N_i+N_gamma+N_delta+j,), False)

    # case without site uncertainty term
    if use_uncert == '2eps':
        for j in range(N_j):
            pars.set_vary((2*N_i+N_gamma+N_delta+2*N_j+N_i+1+j,), False)

    # case without uncertainty terms
    if use_uncert == 'noeps':
        for j in range(N_j+N_i+1):
            pars.set_vary((2*N_i+N_gamma+N_delta+2*N_j+j,), False)

    # if a single station (PURA, class A) should be used as site reference
    if ref_ampl == 'PURA':
        pars.set_vary((2*N_i+N_gamma+N_delta+stas_dict['PURA'],), False)

    return pars


def fas_invert(pkldb=cfg.REFS['REFERENCE_DB'], hcomponent='max', 
               model_name='malagnini_Q0', method='SLSQP', jac='exact',
               use_uncert='3eps', ref_ampl='meanA', subtract_noise=True,
               weights=None, bounds='NE_Italy', save=False):

    DEFAULT_OPTIONS = {'disp': True, 'maxiter': 10000}
    # -----------------------------
    # OPEN INPUT FILE for signal db 
    pklpath = cfg.REFS['PKLSPT_FOLDER'] + '/' + pkldb
    stat_p_s, ml_s, fmin_s, fmax_s, hcompout_s = \
        utils.read_input_component(pklpath, hcomponent, 'signal')

    N_i = stat_p_s.N_ev
    N_j = stat_p_s.N_sta
    N_k = stat_p_s.N_freq
    stas_dict = dict(zip(stat_p_s.stas, list(range(N_j))))

    mag = utils.ml_to_mag_munafo(ml_s)
    m0_calc = utils.mag_to_m0_HK(mag)
    alpha_calc = np.log(m0_calc) 

    # ----------------------------------------------
    # BUILD DATASETS AND PARAMETERS
    # ----------------------------------------------
    # load reference synthetic parameters and build static parameters object
    dirpath = cfg.REFS['SYNTHPARS_FOLDER']

    with open(dirpath + '/' + model_name + '.pkl', 'rb') as f:
        p_model = pickle.load(f)
    N_gamma = p_model.gamma.size
    N_delta = p_model.delta.size

    # naming choice determines characteristics of the model to be used
    run_name = utils.set_run_label(model_name, use_uncert, weights, bounds,
                                   subtract_noise, ref_ampl, method,
                                   hcomponent)

    # fixing pars according to chosen parameter hypothesis (=model)
    p_model = fix_params(p_model, model_name, use_uncert, ref_ampl, stas_dict)

    stat_p = myC.StaticParams(M=stat_p_s.M, F=stat_p_s.F, R=stat_p_s.R, 
                              N_ev=N_i, N_sta=N_j, N_freq=N_k)
    pars_dict = utils.create_pars_dict(p_model, stat_p)

    # -----------------------------------------------
    # update static parameters object - subtract NOISE from data if requested
    if subtract_noise is True:
        Z_true = np.full((N_i*N_j*N_k), 1.)
        stat_p_n, ml_n, fmin_n, fmax_n, hcompout_n = \
            utils.read_input_component(pklpath, hcomponent, 'noise')
        Z_diff = stat_p_s.data - stat_p_n.data
        Z_true[np.where(Z_diff > 0)] = Z_diff[np.where(Z_diff > 0)]
    else:
        Z_true = stat_p_s.data

    stat_p.set_data(Z_true)
    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)

    if bounds == 'NE_Italy':
        inv_bounds = bounds_NE_Italy(p_model, alpha_calc, m0_calc)
    else:
        inv_bounds = bounds

    if bounds is not None:
        if len(bounds) != len(p_model.pars):
            utils.MinimizeException("Bounds must have same lenght as "
                                    "parameter object")

    # -----------------------------------------------
    # set initial parameters from reference
    p_input = deepcopy(p_model)

    # -----------------------------------------------
    # SLSQP REGRESSION - WITH BOUNDS

    if jac == 'exact':
        print('SLSQP, bounded, exact jac')
        start_time = timeit.default_timer()
        if (ref_ampl == 'PURA') or (model_name.startswith('fixA')):
            fmin, p_out = minimizer_wrapper(0, handle_func, p_input, stat_p, 
                                            method='SLSQP',
                                            jac=handle_func_deriv,
                                            options=DEFAULT_OPTIONS, 
                                            bounds=inv_bounds, meanavg=False) 
        else:
            fmin, p_out = minimizer_wrapper(0, handle_func, p_input, stat_p, 
                                            method='SLSQP',
                                            jac=handle_func_deriv,
                                            options=DEFAULT_OPTIONS, 
                                            bounds=inv_bounds)
        elapsed = timeit.default_timer() - start_time
        print('Elapsed time: ', elapsed)
        print(handle_func(p_out, *fcn_args_ee))
        print('------------------')
        if 1: 
            utils.writelog(cfg.REFS['LOG_FOLDER'] + '/' + run_name + '.txt', 
                           fmin, pars_dict, p_input, p_model, p_out)

    else:
        print('SLSQP, bounded, 3-point approx jac')
        start_time = timeit.default_timer()
        fmin, p_out = minimizer_wrapper(0, handle_func, p_input, stat_p, 
                                        method='SLSQP', jac='3-point',
                                        options=DEFAULT_OPTIONS,
                                        bounds=inv_bounds)
        elapsed = timeit.default_timer() - start_time
        print('Elapsed time: ', elapsed)
        print('Cost function of inversion result: ',
              handle_func(p_out, *fcn_args_ee))
        print('------------------')
        if 1:
            utils.writelog(cfg.REFS['LOG_FOLDER'] + '/' + run_name 
                           + '_nojac.txt', fmin, pars_dict, p_input, p_model, 
                           p_out)

    if save is True:
        outpath = cfg.REFS['HOMEPATH'] + '/' + run_name 
        if os.path.exists(outpath) == 0:
            os.makedirs(outpath)
        with open(outpath + '/' + run_name + '.pkl', 'wb') as fout:
            pickle.dump(p_out, fout, protocol=pickle.HIGHEST_PROTOCOL)

    return fmin, p_out
