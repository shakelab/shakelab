import numpy as np


class Params(object):
    """
    Parameters class.
    The class features are:
      alpha <--> ln(M_0) [cm^2 * kg^2 / s^2]
      beta  <--> f_c [Hz]
      gamma <--> gamma indexes of geom. spread. function [adimensional]
      delta <--> Q_0 [adimensional] (if Q_alpha is not used)
      site_fi <--> ln(A) [adimensional]
      site_k  <--> k attenuation factor [s]
      eps_source <--> uncertainty on source contribution
      eps_path <--> uncertainty on path contribution
      eps_site <--> uncertainty on site contribution
      vary <--> flag to fix/release each parameter
      min  <--> lower bound for parameter space
      max  <--> upper bound for parameter space
    """    
    empty = np.array([])

    def __init__(self, alpha=empty, beta=empty, gamma=empty, delta=empty, 
                 site_fi=empty, site_k=empty, eps_source=empty, eps_path=empty, 
                 eps_site=empty):

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.site_fi = site_fi
        self.site_k = site_k
        self.eps_source = eps_source
        self.eps_path = eps_path
        self.eps_site = eps_site
        self.update_p()
        self.vary_all(True)

    def update_p(self):
        self.pars = np.concatenate((self.alpha, self.beta, self.gamma,
                                     self.delta, self.site_fi, self.site_k,
                                     self.eps_source, self.eps_path,
                                     self.eps_site))
        
    def set_alpha(self, a):
        self.alpha = a
        self.update_p()

    def set_beta(self, b):
        self.beta = b
        self.update_p()
        
    def set_gamma(self, c):
        self.gamma = c
        self.update_p()

    def set_delta(self, d):
        self.delta = d
        self.update_p()

    def set_site_fi(self, sf):
        self.site_fi = sf
        self.update_p()

    def set_site_k(self, sk):
        self.site_k = sk
        self.update_p()

    def set_eps_source(self, eso):
        self.eps_source = eso
        self.update_p()

    def set_eps_path(self, epa):
        self.eps_path = epa
        self.update_p()

    def set_eps_site(self, esi):
        self.eps_site = esi
        self.update_p()

    def set_all(self, a, b, c, d, sf, sk, eso, epa, esi):
        self.alpha = a
        self.beta = b
        self.gamma = c
        self.delta = d
        self.site_fi = sf
        self.site_k = sk
        self.eps_source = eso
        self.eps_path = epa
        self.eps_site = esi
        self.update_p()

    def vary_all(self, value):
        self.vary = np.full(len(self.pars), value)
        if value is True:
            self.set_min_all(-1*np.inf)
            self.set_max_all(1*np.inf)

        if value is False:
            for i in range(len(self.pars)):
                self.min[i] = self.pars[i]-1.e-10
                self.max[i] = self.pars[i]+1.e-10

    # NB 'index' must be a list of indexes to change
    def set_vary(self, index, value):
        for i in range(len(index)):
            self.vary[index[i]] = value
            if value is True:
                self.set_min((index[i],), -1*np.inf)
                self.set_max((index[i],), 1*np.inf)

            if value is False:
                self.set_min((index[i],), (self.pars[index[i]]-1.e-10))
                self.set_max((index[i],), (self.pars[index[i]]+1.e-10))

    def set_min(self, index, value):
        if isinstance(value, np.float64) or np.isinf(value):
            for i in range(len(index)):
                self.min[index[i]] = value
        else:
            for i in range(len(index)):
                self.min[index[i]] = value[i]

    def set_max(self, index, value):
        if isinstance(value, np.float64) or np.isinf(value):
            for i in range(len(index)):
                self.max[index[i]] = value
        else:
            for i in range(len(index)):
                self.max[index[i]] = value[i]

    def set_min_all(self, value):
        self.min = np.full(len(self.pars), value)

    def set_max_all(self, value):
        self.max = np.full(len(self.pars), value)


class StaticParams(object):
    """
    Static parameters class.
    The class features are:
      N_ev <--> [scalar] number of events
      N_sta <--> [scalar] number of stations
      N_freq <--> [scalar] number of frequency points (in each spectrum)
      M <--> [array, N_ev*N_sta*N_freq] boolean values to flag
                                        datapoints used for calculations
      F <--> [array, N_ev*N_sta*N_freq] frequency points
      R <--> [array, N_ev*N_sta*N_freq] hypocentral distances
      data <--> [array, N_ev*N_sta*N_freq] FAS data
      weights <--> [array, N_ev*N_sta*N_freq] optional weights applied
                                              to data
      scale  <--> [scalar] optional scaling for cost function
                           calculation and for forward modelling
      orids <--> [array, N_ev] list of event origin IDs
      stas <--> [array, N_sta] list of stations
      freqs <--> [array, N_freqs] list of frequency points in each
                                  spectrum
    """    

    def __init__(self, M=None, F=None, R=None, N_ev=None, N_sta=None, 
                 N_freq=None, data=None, weights=None, scale=None, orids=None,
                 stas=None, freqs=None):

        self.M = M
        self.F = F
        self.R = R
        self.N_ev = N_ev
        self.N_sta = N_sta
        self.N_freq = N_freq
        self.data = data
        self.weights = weights
        self.scale = scale
        self.orids = orids
        self.stas = stas
        self.freqs = freqs
        
    def set_M(self, M):
        self.M = M

    def set_F(self, F):
        self.F = F
        
    def set_R(self, R):
        self.R = R

    def set_N_ev(self, N_ev):
        self.N_ev = N_ev

    def set_N_sta(self, N_sta):
        self.N_sta = N_sta

    def set_N_freq(self, N_freq):
        self.N_freq = N_freq

    def set_data(self, data):
        self.data = data

    def set_scale(self, scale):
        self.scale = scale

    def set_weights(self, weights):
        self.weights = weights

    def set_orids(self, orids):
        self.orids = orids

    def set_stas(self, stas):
        self.stas = stas

    def set_freqs(self, freqs):
        self.freqs = freqs

    def set_all(self, M=None, F=None, R=None, N_ev=None, N_sta=None, 
                N_freq=None, data=None, weights=None, scale=None, orids=None,
                stas=None, freqs=None):
        self.M = M
        self.F = F
        self.R = R
        self.N_ev = N_ev
        self.N_sta = N_sta
        self.N_freq = N_freq
        self.data = data
        self.weights = weights
        self.scale = scale
        self.orids = orids
        self.stas = stas
        self.freqs = freqs


class Spectrum(object):
    """
    Spectrum class (for individual FAS object) containing also info on the event originating it.
    The class features are:
      freq <--> [array, N_freq] frequency points
      amp <--> [array, N_freq] FAS amplitude values
      sm_amp <--> [array, N_freq] smoothed FAS amplitude values
      fbool <--> [array, N_freq] boolean values to flag datapoints used
                                 for calculations
      delta <--> [scalar] source-station angular distance
      meta <--> dictionary containing orid, sta, chan, hypdist
      orid <--> [scalar] origin identifier (event number)
      sta <--> [scalar] station identifier (station name)
      chan <--> [scalar] data channel
      hypdist <--> [scalar] hypocentral distance in km
      ml  <--> [scalar] local magnitude of the recorded event
      freq_lims  <--> [array, 2] lower and upper frequencies used for 
                                 calculations
    """

    freq = np.array([])
    amp = np.array([])
    sm_amp = np.array([])
    fbool = np.array([])
    delta = []
    meta = {}
    orid = " "
    sta = " "
    chan = " "
    hypdist = " "
    ml = " "
    freq_lims = np.array([0., 0.])

    def __init__(self, orid, sta, chan, hypdist, ml, f1, f2):
        self.orid = orid
        self.sta = sta
        self.chan = chan
        self.hypdist = hypdist
        self.ml = ml
        self.freq_lims[0] = f1
        self.freq_lims[1] = f2
        self.meta = self.__set_trace_meta()

    def __set_trace_meta(self):
        nm = {}
        nm.update({self.orid: 'orid'})
        nm.update({self.sta: 'sta'})
        nm.update({'chan': self.chan})
        nm.update({'hypdist': self.hypdist})
        return nm

    def set_trace_freq(self, f):
        self.freq = f

    def set_trace_amp(self, a):
        self.amp = a

    def set_trace_smamp(self, sa):
        self.sm_amp = sa

    def set_trace_fbool(self, b):
        self.fbool = b

    def set_trace_delta(self, d):
        self.delta = d

    def set_trace_flim(self, f1, f2):
        self.freq_lims = [f1, f2]

    def __str__(self):
        string = u'[<Spectrum> orid:%s sta:%s chan:%s hypdist:%s ml:%s ' \
                 u'f1:%s f2:%s \n freq:%s \n amp:%s \n sm_amp:%s \n fbool:%s' \
                 u' \n delta:%s]' % (self.orid, self.sta, self.chan,
                                     self.hypdist, self.ml, self.freq_lims[0],
                                     self.freq_lims[1], self.freq, self.amp,
                                     self.sm_amp, self.fbool, self.delta)
        return string
