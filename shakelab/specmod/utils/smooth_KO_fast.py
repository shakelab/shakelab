import numpy as np

def _smooth_konno_ohmachi_fast(frequencies, spectrum, fcs, bandwidth=40.):
    n = 3
    upper_limit = np.power(10, +n/bandwidth)
    lower_limit = np.power(10, -n/bandwidth)
    nrows = 1
    ncols = fcs.size
    smoothed_spectrum = np.empty((nrows, ncols))
    for fc_index, fc in enumerate(fcs):
        if fc < 1E-6:
            smoothed_spectrum[:,fc_index] = 0
            continue
        sumproduct = np.zeros(nrows)
        sumwindow = 0
        for f_index, f in enumerate(frequencies):
            f_on_fc = f/fc
            if (f < 1E-6) or (f_on_fc > upper_limit) or (f_on_fc < lower_limit):
                continue
            elif np.abs(f - fc) < 1e-6:
                window = 1.
            else:
                window = bandwidth * np.log10(f_on_fc)
                window = np.sin(window) / window
                window *= window
                window *= window
            sumproduct += window*spectrum[f_index]
            sumwindow += window
        if sumwindow > 0:
            smoothed_spectrum[:,fc_index] = sumproduct / sumwindow
        else:
            smoothed_spectrum[:,fc_index] = 0
    return smoothed_spectrum[0]