import shutil
import os
import numpy as np
import pickle
from matplotlib import pyplot as plt
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.utils as utils
from spectral_modelling.model_and_invert.fas_log import ln_fas


def avg_mean_residual(stat_p, Z_calc):
    """
    Calculates freq.-dependent site amplification as geometric average of the
    inversion residuals (GAVG) and the related geometric standard deviation 
    (GSD)
    """

    N_i = stat_p.N_ev
    N_j = stat_p.N_sta
    N_k = stat_p.N_freq
    GAVG = np.ones((N_j * N_k))
    GSD = np.ones((N_j * N_k))

    for s in range(N_j):
        i = 0
        mask_matrix = np.zeros((N_i, N_k))
        aj_matrix = np.ones((N_i, N_k))

        for e in range(N_i):
            F = np.zeros(N_k)
            Z = np.zeros(N_k)
            M = np.zeros(N_k)
            Z_c = np.zeros(N_k)
            for k in range(N_k):
                F[k] = stat_p.F[k + N_k * s + N_j * N_k * e]
                Z[k] = stat_p.data[k + N_k * s + N_j * N_k * e]
                M[k] = stat_p.M[k + N_k * s + N_j * N_k * e]
                Z_c[k] = Z_calc[k + N_k * s + N_j * N_k * e]
            index = np.where(M == 1)
            if len(index) > 0:
                aj_matrix[i][index] = (Z / Z_c)[index]
                mask_matrix[i][:] = M
                i += 1

        mask = np.sum(mask_matrix, axis=0)
        index5 = np.where(mask >= 5)
        mask_5 = mask[index5]

        gmean_5 = np.ones_like(mask_5)
        gsd_5 = np.zeros_like(mask_5)

        aj_matrix_5 = aj_matrix[:, index5[0]]
        for k in range(len(mask_5)):
            aa = aj_matrix_5[:, k]
            bb = aa[np.where(aa != 1.)]
            gmean_5[k] = np.power(np.prod(bb, axis=0), 1. / mask_5[k])

        freqs_5 = stat_p.freqs[index5]

        # NOTE maybe join with previous cycle
        for k in range(len(mask_5)):
            aa = aj_matrix_5[:, k]
            bb = aa[np.where(aa != 1.)]
            gsd_5[k] = np.sum((np.log(bb / gmean_5[k])) ** 2, axis=0) / \
                       mask_5[k]
            gsd_5[k] = np.exp(np.sqrt(gsd_5[k]))

        for k in range(N_k):
            f = stat_p.freqs[k]
            if f in freqs_5:
                ind_f = np.where(freqs_5 == f)
                GAVG[k + N_k * s] = gmean_5[ind_f][0]
                GSD[k + N_k * s] = gsd_5[ind_f][0]

        return GAVG, GSD


def fdsiteamp(run_name, pkldb=cfg.REFS['REFERENCE_DB'], hcomponent='max',
              plot=True):
    """
    This feature calculates and saves the frequency-dependent site 
    amplification functions a(f) obtained from residuals.
    Optionally, it creates plots for a(f) and A*a(f) functions.
    """

    # OPEN INPUT FILEs
    pklpath = cfg.REFS['PKLSPT_FOLDER'] + '/' + pkldb
    stat_p_s, ml_s, fmin_s, fmax_s, compout_s = \
        utils.read_input_component(pklpath, hcomponent)
    stat_p_n, ml_n, fmin_n, fmax_n, compout_n = \
        utils.read_input_component(pklpath, hcomponent, 'noise')

    N_i = stat_p_s.N_ev
    N_j = stat_p_s.N_sta
    N_k = stat_p_s.N_freq
    orids_dict = dict(zip(stat_p_s.orids, list(range(N_i))))
    stas_dict = dict(zip(stat_p_s.stas, list(range(N_j))))
    reverse_stas_dict = dict(zip(list(range(N_j)), stat_p_s.stas))

    model_name = utils.read_run_label(run_name)[1]
    dirpath = cfg.REFS['SYNTHPARS_FOLDER']
    with open(dirpath + '/' + model_name + '.pkl', 'rb') as f:
        p_model = pickle.load(f)

    N_gamma = p_model.gamma.size
    N_delta = p_model.delta.size

    stat_p = stat_p_s

    inpath = cfg.REFS['HOMEPATH'] + '/' + run_name 
    with open(inpath + '/' + run_name + '.pkl', 'rb') as f:
        Params_slsqp = pickle.load(f)

    # extract useful parameters from Params_slsqp
    site_fi_inv = Params_slsqp.site_fi

    # build the Params object equivalent to Params_slsqp, but without eps,
    # to be used to calculate the forward model
    p_strip = utils.strip_pars(Params_slsqp, stat_p)

    ee = np.array([])
    fcn_args_ee = (stat_p, ee, ee, ee, N_gamma, N_delta)
    Z_calc = ln_fas(p_strip, fcn_args_ee)
    Z_calc = np.exp(Z_calc)

    Z_true = np.full((N_i * N_j * N_k), 1.)
    Z_diff = stat_p_s.data - stat_p_n.data
    Z_true[np.where(Z_diff > 0)] = Z_diff[np.where(Z_diff > 0)]

    stat_p.set_data(Z_true)

    homepath = cfg.REFS['HOMEPATH'] + '/' + run_name 
    if os.path.exists(homepath) == 0:
        os.makedirs(homepath)

    GAVG, GSD = avg_mean_residual(stat_p, Z_calc)

    #################
    # save mean residuals (f, a(f), associated sigma) for each station
    dirpath = homepath + '/meanresiduals'
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)
    os.makedirs(dirpath)

    for s in range(len(stas_dict)):
        filename = dirpath + '/' + str(reverse_stas_dict[s]) + '_siteamp.txt'
        with open(filename, 'w') as fout:
            for k in range(N_k):
                fout.write(str(stat_p.freqs[k]) + ' ' + str(
                    GAVG[k + N_k * s]) + ' ' + str(GSD[k + N_k * s]) + '\n')

    #################
    # plot mean residuals (with variation) for each station
    # and total site amplification (site_fi * meanresidual)
    if plot is True:

        dirpath = homepath + '/plots/meanresiduals'
        if os.path.exists(dirpath):
            shutil.rmtree(dirpath)
        os.makedirs(dirpath)

        for s in range(len(stas_dict)):

            fig, ax = plt.subplots(1, 2, figsize=(6, 3))
            fig.suptitle(
                'sta = ' + str(reverse_stas_dict[s]) + ' , soil class = '
                + str(_c.soil_dict[reverse_stas_dict[s]]) + ', '
                + hcomponent + ', ' + run_name, fontsize=8)

            for e in range(len(orids_dict)):
                F_e = np.zeros(N_k)
                Z_e = np.zeros(N_k)
                M_e = np.zeros(N_k)
                Z_calc_e = np.zeros(N_k)
                for k in range(N_k):
                    F_e[k] = stat_p.F[k + N_k * s + N_j * N_k * e]
                    Z_e[k] = Z_true[k + N_k * s + N_j * N_k * e]
                    M_e[k] = stat_p.M[k + N_k * s + N_j * N_k * e]
                    Z_calc_e[k] = Z_calc[k + N_k * s + N_j * N_k * e]
                index = np.where(M_e == 1)
                F_1 = F_e[index]
                Z_1 = Z_e[index]
                Z_CALC_1 = Z_calc_e[index]
                if len(F_1) > 0:
                    ax[0].loglog(F_1, (Z_1 / Z_CALC_1), c='grey', lw=0.5)
                    ax[1].loglog(F_1, (Z_1 / Z_CALC_1) * np.exp(site_fi_inv[s]),
                                 c='grey', lw=0.5)

            gavg_s = GAVG[(N_k * s): (N_k * s + N_k)]
            gsd_s = GSD[(N_k * s): (N_k * s + N_k)]
            index = np.where(gavg_s != 1.)
            ax[0].loglog(stat_p.freqs[index], gavg_s[index], c='blue', zorder=2)
            ax[0].loglog(stat_p.freqs[index], (gavg_s * gsd_s)[index], lw=0.8,
                         c='blue', zorder=2)
            ax[0].loglog(stat_p.freqs[index], (gavg_s / gsd_s)[index], lw=0.8,
                         c='blue', zorder=2)
            ax[1].loglog(stat_p.freqs[index],
                         gavg_s[index] * np.exp(site_fi_inv[s]), c='blue',
                         zorder=2)

            for i in range(2):
                ax[i].axhline(y=1., zorder=1, c='k', lw=0.8)
                ax[i].axhline(y=2.5, zorder=1, c='k', lw=0.8, ls='dotted')
                ax[i].axhline(y=0.4, zorder=1, c='k', lw=0.8, ls='dotted')
                ax[i].set_xlabel('frequency [Hz]', fontsize=7)
                ax[i].tick_params(axis='both', which='both', labelsize=7)
                ax[i].set_ylim([0.05, 40])

            ax[0].set_ylabel('residuals (FAS_real - FAS_modelled)', fontsize=7)
            ax[0].set_title('freq. dependent amplification', fontsize=7)
            ax[1].set_ylabel('residuals (FAS_real - FAS_modelled) * A',
                             fontsize=7)
            ax[1].set_title('scaled freq. dependent amplification', fontsize=7)

            plt.subplots_adjust(wspace=0.35, top=0.85, right=0.97, left=0.1,
                                bottom=0.15)
            plt.savefig(dirpath + '/' + str(reverse_stas_dict[s]) + '.png')
            plt.close(fig)
