import numpy as np
import pickle
from shakelab.libutils.geodetic import tunnel_distance_sphere
from spectral_modelling.utils.myClass import Spectrum
from spectral_modelling.utils.smooth_KO_fast import _smooth_konno_ohmachi_fast
import scipy.fftpack as fftpack
import spectral_modelling.utils.constants as _c
import spectral_modelling.utils.config as cfg
import spectral_modelling.utils.utils as utils


inpath = cfg.REFS['INPUT_FOLDER']
originmeta = cfg.REFS['ORIGINMETA']
stationmeta = cfg.REFS['STATIONMETA']
originlist = cfg.REFS['ORIGINDB']
stationlist = cfg.REFS['STATIONDB']

#TBD: add cutting the waveform in a desired interval accoridng to P or S picks

def create_spectralpkl():
    """
    Read processed waveforms for selected events and stations, calculate
    their smoothed fft spectra, and store inofrmation separately for each of the horizontal components in form of pickle files containing a list of 
    Spectrum() objects.

    The waveforms are read from a pickle containing a shakelab 
    StreamCollection() object.
    """

    # parse metadata for all events
    evmetadb = np.genfromtxt(inpath + '/' + originmeta, names=True)
    allids = evmetadb['orid']
    allids = [int(i) for i in allids]
    mags = dict(zip(allids, evmetadb['ml']))
    elats = dict(zip(allids, evmetadb['evlat']))
    elons = dict(zip(allids, evmetadb['evlon']))
    edepths = dict(zip(allids, evmetadb['hypdep']))

    # parse metadata for all stations
    stasmetadb = np.genfromtxt(inpath + '/' + stationmeta, names=True, 
                               dtype=None, encoding='UTF-8')
    allstas = stasmetadb['stacode']
    slats = dict(zip(allstas, stasmetadb['slat']))
    slons = dict(zip(allstas, stasmetadb['slon']))
    selev = dict(zip(allstas, stasmetadb['selev']))

    # parse id of selected events
    evdb = np.genfromtxt(inpath + '/' + originlist, names=True)
    orids = [int(i) for i in evdb['orid']]

    # parse id of selected stations
    stasdb = np.genfromtxt(inpath + '/' + stationlist, names=True, dtype=None, encoding='UTF-8')
    stas = stasdb['stacode']

    outpath = cfg.REFS['PKLSPT_FOLDER'] + '/' + cfg.REFS['REFERENCE_DB']
    out_EW = open(outpath + 'EW.pkl', 'wb')
    out_NS = open(outpath + 'NS.pkl', 'wb')

    # Loop over selected orids
    for orid in orids:
        ml = mags[orid]
        elat = elats[orid]
        elon = elons[orid]
        edepth = edepths[orid]
        pklname = cfg.REFS['PKLMSEED_FOLDER'] + '/' + str(orid) + '.pkl'
        with open(pklname, 'rb') as f:
            evstreams = pickle.load(f)

        # Loop over selected stations
        sids = evstreams.sid
        sids = [k for k in sids if (k.split('.')[1] in stas) and 
                (k.endswith(tuple(['E','N'])))]
        for sid in sids:
            tr = evstreams[sid].get()
            sta = utils.nslc(sid)['sta']
            chan = utils.nslc(sid)['chan']
            orientation = chan[-1]

            hypdist = tunnel_distance_sphere(elat, elon, edepth, slats[sta], 
                                            slons[sta], selev[sta])
            #per ora fisse
            f1 = float('0.5')
            f2 = float('25.0')
            w = Spectrum(orid, sta, chan, hypdist, ml, f1, f2)
            w.set_trace_flim(f1, f2)

            npts = tr.head.nsamp
            dt = tr.head.delta

            # Calculate amplitude fft spectra
            nfrqs = int(npts/2)+1 if (npts % 2) == 0 else int((npts+1)/2)
            frq = np.abs(np.fft.fftfreq(npts, dt))[:nfrqs]
            hw = np.array([(1 - np.cos(2 * _c.PI * j / npts))/2. 
                           for j in range(npts)])
            tr_fft = np.sqrt(2*((np.abs(fftpack.fft(tr.data)[:nfrqs]))**2)
                             / (np.sum(hw**2)) / dt)
            # Calculate smoothed fft
            fcs = cfg.MODELS["FREQS"]
            fcs = np.array(fcs, dtype=np.double)
            sm_fft = _smooth_konno_ohmachi_fast(frq, tr_fft, fcs, bandwidth=40.)

            w.set_trace_freq(fcs)
            w.set_trace_smamp(sm_fft)

            fbool = np.zeros_like(fcs)
            fbool[np.where(((fcs-f1) >= 1.e-04) & ((fcs-f2) <= 1.e-04))] = 1
            w.set_trace_fbool(fbool)

            match orientation:
                case 'E':
                    pickle.dump(w, out_EW, protocol=pickle.HIGHEST_PROTOCOL)
                case 'N':
                    pickle.dump(w, out_NS, protocol=pickle.HIGHEST_PROTOCOL)

    out_EW.close()
    out_NS.close()

