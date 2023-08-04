"""
Constants needed in FAS calculation
"""

import math

PI = math.pi
EE = math.e

RHO = 2.8  # g/cm^3
V_S = 350000  # cm/s

R0 = 100000.  # cm
Rref = 7000000.  # cm

R1 = 4000000.  # cm
R2 = 5000000.  # cm
R3 = 6000000.  # cm
R4 = 8000000.  # cm
R5 = 10000000.  # cm

R_THETA = 0.55
FF = 2.
CSI = 1. / math.sqrt(2.)
LN_MCOST = math.log((R_THETA * FF * CSI) / (4. * PI * RHO * (V_S ** 3) * R0))
MCOST = (R_THETA * FF * CSI) / (4. * PI * RHO * (V_S ** 3) * R0)

QCOST = PI / V_S

AVG_STDROP = 7.3e+06  # g/(cm*s^2)

soil_dict = {'AVS': 'A', 'CARC': 'C', 'CESC': 'B', 'CMO': 'A', 'DST2': 'A',
             'FDS': 'A', 'FLP': 'C', 'GEDE': 'B', 'GEPF': 'A', 'GESC': 'C',
             'MASA': 'A', 'MOGG': 'A', 'PRAD': 'A', 'RST': 'A', 'SPP': 'C',
             'AUP': 'A', 'CHF': 'A', 'PAUL': 'A', 'POLC': 'B', 'PURA': 'A',
             'STOL': 'C', 'TLM2': 'B', 'GORI': 'B', 'DANT': 'A'}
