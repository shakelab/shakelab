import numpy as np
import scipy.optimize as opt
from copy import deepcopy
import spectral_modelling.utils.utils as utils
import spectral_modelling.utils.config as cfg

# -----------------------------------------------------------------------------
# Wrapper to call scipy.optimize.minimize without having to change the
# input parameter to an array, and to refurbish the output as a Param()
# item
# -----------------------------------------------------------------------------


def minimizer_tool(fcn, params_in, fcn_args, bounds=None, method='CG',
                   jac=None, tol=1.e-08, options=None, diff_step=None,
                   jac_sparsity=None, meanavg=True):
    """Perform a fit of a set of parameters by minimizing an objective
    (cost) function using one of the several available methods in scipy.

    The minimize function takes an objective function to be minimized,
    a dictionary (:class:`~myClass.Params`) containing the model
    parameters, and several optional arguments.
    """

    stat_params, index_v, p_fixed, index_f, N_gamma, N_delta = fcn_args
    utils.checkext_statparams(stat_params)
    utils.checkext_params(params_in)
    pars = params_in.pars

    N_ev = stat_params.N_ev
    N_sta = stat_params.N_sta

    # define freq.-independent amplification relative to the network
    # average
    # NB1 now it uses the mean avg on rock (class A), as defined in
    #     utils.py
    # NB2 the corresponding jacobian can also be calculated and added to
    #     the site_avg dictionary as 'jac' (optional)

    def constraint_A(amps_in, N_ev, N_sta, N_gamma, N_delta, index_f, p_fixed):
        amps_input = utils.insert_elem(amps_in, index_f, p_fixed)      
        flatindex = utils.stas_flat()
        amps_sum = 0.
        for i in flatindex:
            amps_sum = amps_sum + amps_input[2*N_ev+N_gamma+N_delta+i]
        return amps_sum

    site_avg = ({'type': 'eq',
                 'fun': constraint_A,
                 'args': (N_ev, N_sta, N_gamma, N_delta, index_f, p_fixed)
                 })

    # -----------------------------------------
    # perform minimization

    if method in cfg.MINIMIZE_METHODS_BOUNDS:
        mybounds = opt.Bounds(np.array(bounds[0]), np.array(bounds[1]))
        if method == 'SLSQP':
            if meanavg is True:
                result = opt.minimize(fcn, pars, args=fcn_args,
                                      bounds=mybounds, method=method, jac=jac,
                                      constraints=site_avg, tol=tol,
                                      options=options)
            else:
                result = opt.minimize(fcn, pars, args=fcn_args,
                                      bounds=mybounds, method=method, jac=jac,
                                      tol=tol, options=options)

        else:
            result = opt.minimize(fcn, pars, args=fcn_args, bounds=mybounds, 
                                  method=method, jac=jac, tol=tol,
                                  options=options)

    elif method in cfg.MINIMIZE_METHODS_NOBOUNDS:
        result = opt.minimize(fcn, pars, args=fcn_args, method=method, jac=jac, 
                              tol=tol, options=options)            
       
    elif method == "least_squares":
        if jac is None:
            jac = '2-point'           
            
        if diff_step is not None:
            result = opt.least_squares(fcn, pars, args=fcn_args, bounds=bounds, 
                                       jac=jac, ftol=tol, diff_step=diff_step, 
                                       verbose=2)

        if jac_sparsity is not None:
            result = opt.least_squares(fcn, pars, args=fcn_args, bounds=bounds, 
                                       jac_sparsity=jac_sparsity, ftol=tol, 
                                       verbose=2)

        else:
            result = opt.least_squares(fcn, pars, args=fcn_args, bounds=bounds, 
                                       jac=jac, ftol=tol, verbose=2)

    else:
        raise utils.MinimizeException("Method not recognized")
    
    # -----------------------------------------
    # build output Params object

    xout = result.x
    xout_all = utils.insert_elem(xout, index_f, p_fixed)

    params_out = utils.array_to_Param(xout_all, stat_params, N_gamma, N_delta)

    return result, params_out
       

def minimizer_wrapper(index, fcn, p_input, stat_p, method='CG', jac=None, 
                      tol=1.e-08, options=None, diff_step=None, 
                      jac_sparsity=None, bounds=None, meanavg=True):
    """Actual wrapper; uses bounds (when applicable) to keep fixed
       parameters actually locked.

       NB bounds must be either a list or an array with two components,
       the first for lower bound values and the second for upper bound
       values
    """

    utils.checkext(p_input, stat_p)

    p_in, index_v, p_fixed, index_f = utils.create_p_in(p_input, stat_p)

    N_gamma = p_input.gamma.size
    N_delta = p_input.delta.size
    fcn_args = (stat_p, index_v, p_fixed, index_f, N_gamma, N_delta)

    # manually 'lock' fixed values by using equal bounds
    start_vals = deepcopy(p_in)

    if bounds is None:
        lower_bounds, upper_bounds = [], []
        for i in range(len(p_in.pars)):
            par = p_in.pars[i]
            lower_bounds.append(-1*np.inf)
            upper_bounds.append(1*np.inf)

    else:
        lower_bounds, upper_bounds = [], []
        lbounds_v = list(np.array(bounds[0])[index_v])
        ubounds_v = list(np.array(bounds[1])[index_v])
        for i in range(len(p_in.pars)):
            lower_bounds.append(lbounds_v[i])
            upper_bounds.append(ubounds_v[i])

    res, p_out = minimizer_tool(fcn, start_vals, fcn_args, 
                                bounds=(lower_bounds, upper_bounds), 
                                method=method, jac=jac, tol=tol, 
                                options=options, diff_step=diff_step, 
                                jac_sparsity=jac_sparsity, meanavg=meanavg)

    p_out.vary_all(True)
    p_out.set_vary(index_f, False)

    return res, p_out
