# Configuration File 
# -----
# SRCMODEL = model for source contribution
# MOTION = motion component (acc/vel/disp)
# FREQS = frequency points at which spectra should be calculated
# -----
# HOMEPATH = home folder path for calculations
# SPTPATH = path to the folder containing input .spt files
# SYNTHPARS_FOLDER = name of the folder storing synthetic (initial) 
#                    parameter sets
# TXTINPUT_FOLDER = name of the folder containing additional input .txt
#                   files
# ORIGIN_DB = name of the txt file containing reference origins 
#             database information (located in TXTINPUT_FOLDER)
# REFERENCE_DB = name of pkl database created by parser_to_cpickle.py
#                and used in utils.py to get the set of stations in
#                EC8-class A (used to set the reference condition on
#                site_fi)
# REFRENCE_DB_VERSION = version of pkl database (useful if more
#                       databases are extracted from same underlying
#                       data)

import numpy as _np

# lists of recognized minimazion methods
MINIMIZE_METHODS_BOUNDS = ['Powell', 
                           'L-BFGS-B', 
                           'TNC', 
                           'trust-constr',
                           'SLSQP']

MINIMIZE_METHODS_NOBOUNDS = ['Nelder-Mead', 
                             'Powell', 
                             'CG', 
                             'BFGS',
                             'Newton-CG', 
                             'L-BFGS-B', 
                             'TNC', 
                             'COBYLA',
                             'SLSQP', 
                             'trust-constr', 
                             'dogleg', 
                             'trust-ncg',
                             'trust-exact', 
                             'trust-krylov']

# modelling choices
MODELS = dict(
    SRCMODEL='BRUNE',
    MOTION='velocity',
    FREQS=_np.logspace(_np.log10(0.5), _np.log10(25), num=30, endpoint=True)
)

# reference naming and path choices
REFS = dict(
    HOMEPATH='/path/to/homepath',
    INPUT_FOLDER='/path/to/inputfiles',
    LOG_FOLDER='/path/to/log',
    SYNTHPARS_FOLDER='/path/to/synth_params',
    PKLMSEED_FOLDER='/path/to/event_mseeds',
    PKLSPT_FOLDER='/path/to/event_spectra',
    REFERENCE_DB='vel_fullwave_',
    ORIGINMETA='origins_coordinates.txt',
    STATIONMETA='stations_coordinates.txt',
    ORIGINDB='filt_origins_list.txt',
    STATIONDB='filt_stations_list.txt'
)

