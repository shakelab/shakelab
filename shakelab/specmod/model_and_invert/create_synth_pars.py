import numpy as np
import pickle
import os
import shutil
import spectral_modelling.utils.utils as utils
import spectral_modelling.utils.myClass as myC
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.constants as _c

def synthetic_p_in(dbname=cfg.REFS['REFERENCE_DB'], hcomponent='max', 
                   model='malagnini_Q0'):
    """
    Create sample input parameter set(s) and store it as pickle object
    """

    # OPEN INPUT FILE and read information
    pklpath = cfg.REFS['PKLSPT_FOLDER'] + '/' + dbname
    stat_p, ml, fmin, fmax, hcompout = utils.read_input_component(pklpath, 
                                                                  hcomponent)

    N_i = stat_p.N_ev
    N_j = stat_p.N_sta

    mag = utils.ml_to_mag_munafo(ml)
    m0_calc = utils.mag_to_m0_HK(mag)
    alpha_calc = np.log(m0_calc)
    beta_calc = utils.m0_to_beta_Brune(m0_calc, _c.AVG_STDROP)

    dirpath = cfg.REFS['SYNTHPARS_FOLDER']
    if os.path.exists(dirpath) == 0:
        os.makedirs(dirpath)

    # build and save a set of model parameters
    # "malagnini_Q0": fc from stress drop; deltas from Malagnini;
    #                 delta[1] fixed to 0
    match model:
        case 'malagnini_Q0':
            alpha = alpha_calc
            beta = beta_calc
            gamma = np.array([-1., -1.6, -1.2, -1.3, -0.5, -0.95, -1.2, -1.8, 
                              -1.2, -0.5], dtype='d')
            delta = np.array([300., 0.], dtype='d')
            site_fi = np.full(N_j, 0.)
            site_k = np.full(N_j, 0.04)
            eps_source = np.zeros(N_i)
            eps_path = np.array([0.])
            eps_site = np.zeros(N_j)
            pars = myC.Params(alpha, beta, gamma, delta, site_fi, site_k, 
                              eps_source, eps_path, eps_site)

            N_gamma = gamma.size
            for i in range(N_gamma):
                pars.set_vary((2 * N_i + i,), False)
            pars.set_vary((2 * N_i + N_gamma + 1,), False)

        case _:
            raise ValueError(f"Model {model} not implemented")

    with open(dirpath + '/' + str(model) + '.pkl', 'wb') as out:
        pickle.dump(pars, out, protocol=pickle.HIGHEST_PROTOCOL)
