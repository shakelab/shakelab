import numpy as np

def Weichert(N, IT, FMAG, pas):
    """
    Weichert method for estimating the Gutenberg-Richter b-value.

    Args:
        N (list): List of earthquake counts per magnitude bin.
        IT (list): List of completeness values for each magnitude bin.
        FMAG (list): List of magnitude values corresponding to the bins.
        pas (float): Increment for magnitude bins.

    Returns:
        tuple: A tuple containing the following values:
            - FMAG (list): Updated magnitude values.
            - mo (float): Value of the minimum magnitude.
            - BETA (float): Estimated b-value.
            - STDV (float): Standard deviation of the b-value.
            - B (float): Estimated b-value in base 10 logarithm.
            - STDB (float): Standard deviation of the b-value in base 10 logarithm.
            - FNGTMO (float): Estimated cumulative number of earthquakes above the completeness magnitude.
            - STDFNGTMO (float): Standard deviation of FNGTMO.

    This function applies the Weichert method to estimate the Gutenberg-Richter b-value,
    which characterizes the distribution of earthquake magnitudes in a seismic region.

    The method iteratively refines the b-value estimate until convergence.

    Note:
    The code assumes 0-based indexing for the bins.

    """
    LOW = 0  # Python uses 0-based indexing
    HIGH = len(N)

    SNM = 0
    NKNOUT = 0
    STMEX = 0
    SUMTEX = 0
    STM2X = 0
    SUMEXP = 0

    BETA = 1.0 / np.log(10)  # initial value

    i = 0
    for k in range(LOW, HIGH):  # 0 till N-1 of bins
        SNM += N[k] * FMAG[k]
        NKNOUT += N[k]
        TJEXP = IT[k] * np.exp(-BETA * FMAG[k])
        TMEXP = TJEXP * FMAG[k]
        SUMEXP += np.exp(-BETA * FMAG[k])
        STMEX += TMEXP
        SUMTEX += TJEXP
        STM2X += FMAG[k] * TMEXP

    V_STMEX = [STMEX]

    DLDB = STMEX / SUMTEX
    D2LDB2 = NKNOUT * (DLDB ** 2 - STM2X / SUMTEX)
    DLDB = DLDB * NKNOUT - SNM
    BETL = BETA

    BETA = BETA - DLDB / D2LDB2
    STDV = np.sqrt(-1.0 / D2LDB2)
    B = BETA / np.log(10)
    STDB = STDV / np.log(10)
    FNGTMO = NKNOUT * SUMEXP / SUMTEX

    while abs(BETA - BETL) >= 0.0001:
        i += 1
        SNM = 0
        NKNOUT = 0
        STMEX = 0
        SUMTEX = 0
        STM2X = 0
        SUMEXP = 0

        for k in range(LOW, HIGH):
            SNM += N[k] * FMAG[k]
            NKNOUT += N[k]
            TJEXP = IT[k] * np.exp(-BETA * FMAG[k])
            TMEXP = TJEXP * FMAG[k]
            SUMEXP += np.exp(-BETA * FMAG[k])
            STMEX += TMEXP
            SUMTEX += TJEXP
            STM2X += FMAG[k] * TMEXP

        DLDB = STMEX / SUMTEX
        D2LDB2 = NKNOUT * (DLDB ** 2 - STM2X / SUMTEX)
        DLDB = DLDB * NKNOUT - SNM
        BETL = BETA
        BETA = BETA - DLDB / D2LDB2
        STDV = np.sqrt(-1.0 / D2LDB2)
        B = BETA / np.log(10)
        STDB = STDV / np.log(10)
        FNGTMO = NKNOUT * SUMEXP / SUMTEX

        V_STMEX.append(STMEX)

    STDFNGTMO = np.sqrt(FNGTMO / NKNOUT)

    mo = FMAG[LOW] - pas / 2

    return FMAG, mo, BETA, STDV, B, STDB, FNGTMO, STDFNGTMO





import numpy as np
from scipy.optimize import minimize

def weichert_objective(BETA, N, IT, FMAG):
    # Objective function for the Weichert method
    TJEXP = IT * np.exp(-BETA * FMAG)
    TMEXP = TJEXP * FMAG
    SUMTEX = np.sum(TJEXP)
    STMEX = np.sum(TMEXP)
    SNM = np.sum(N * FMAG)
    D2LDB2 = np.sum(N * (FMAG - (SNM / SUMTEX)) ** 2) - np.sum(FMAG * N * (FMAG - (SNM / SUMTEX)) / SUMTEX)
    return D2LDB2

def Weichert(N, IT, FMAG, pas):
    """
    Weichert method for estimating the Gutenberg-Richter b-value.

    Args:
        N (numpy.ndarray): Array of earthquake counts per magnitude bin.
        IT (numpy.ndarray): Array of completeness values for each magnitude bin.
        FMAG (numpy.ndarray): Array of magnitude values corresponding to the bins.
        pas (float): Increment for magnitude bins.

    Returns:
        tuple: A tuple containing the following values:
            - FMAG (numpy.ndarray): Updated magnitude values.
            - mo (float): Value of the minimum magnitude.
            - BETA (float): Estimated b-value.
            - STDV (float): Standard deviation of the b-value.
            - B (float): Estimated b-value in base 10 logarithm.
            - STDB (float): Standard deviation of the b-value in base 10 logarithm.
            - FNGTMO (float): Estimated cumulative number of earthquakes above the completeness magnitude.
            - STDFNGTMO (float): Standard deviation of FNGTMO.

    This function applies the Weichert method to estimate the Gutenberg-Richter b-value.

    Note:
    The code assumes 0-based indexing for the bins.

    """
    LOW = 0  # Python uses 0-based indexing
    HIGH = len(N)

    BETA_0 = 1.0 / np.log(10)  # initial value for BETA

    result = minimize(weichert_objective, BETA_0, args=(N, IT, FMAG), method='BFGS', options={'disp': False})
    BETA = result.x[0]

    SNM = np.sum(N * FMAG)
    NKNOUT = np.sum(N)
    TJEXP = IT * np.exp(-BETA * FMAG)
    TMEXP = TJEXP * FMAG
    SUMTEX = np.sum(TJEXP)
    STMEX = np.sum(TMEXP)
    STDV = np.sqrt(1.0 / result.hess_inv[0, 0])
    B = BETA / np.log(10)
    STDB = STDV / np.log(10)
    FNGTMO = NKNOUT * np.sum(np.exp(-BETA * FMAG)) / SUMTEX
    STDFNGTMO = np.sqrt(FNGTMO / NKNOUT)

    mo = FMAG[LOW] - pas / 2

    return FMAG, mo, BETA, STDV, B, STDB, FNGTMO, STDFNGTMO
