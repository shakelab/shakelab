# ****************************************************************************
# Python porting of the original fortran code by F.J. Sanchez Sesma
#
# This file is part of ShakeLab.
# Last modified: 22/7/2022
#
# NOTE: this is a direct porting from the original fortran 77, therefore
# little optimisation is done (e.g. no vectorisation)
# ****************************************************************************

import numpy as np
import cmath as cm


def psvq_soil_response(freq, hl, vp, vs, dn, qp, qs, iwave='sh', iangle=0.):
    """
    Generic interface to PSVQ functions.
    """
    lnum = len(hl)
    fnum = len(freq)
    
    if iwave == 'sh':
        sh_h = np.zeros((lnum, fnum), dtype=complex)

        for i in range(fnum):

            oh = hassh(freq[i], hl, vs, dn, qs, iangle)
            sh_h[:, i] = oh.conj()

        return sh_h

    elif iwave == 'p':
        p_h = np.zeros((lnum, fnum), dtype=complex)
        p_v = np.zeros((lnum, fnum), dtype=complex)

        for i in range(fnum):
            ov, oh = haspsv(freq[i], hl, vp, vs, dn, qp, qs, 1, iangle)
            p_v[:, i] = ov.conj()
            p_h[:, i] = oh.conj()

        return p_v, p_h

    elif iwave == 'sv':
        sv_h = np.zeros((lnum, fnum), dtype=complex)
        sv_v = np.zeros((lnum, fnum), dtype=complex)

        for i in range(fnum):
            ov, oh = haspsv(freq[i], hl, vp, vs, dn, qp, qs, 2, iangle)
            sv_v[:, i] = ov.conj()
            sv_h[:, i] = oh.conj()

        return sv_v, sv_h

    else:
        raise ValueError('Not a valid iwave type')

def hassh(FREQ, H, BETA, RHO, QS, GAMMA=0., POL=0., IDT=-1):
    """
    INPUT
        FREQ : calculation frequency
        H : layer thickness
        BETA : S-wave velocity profile
        RHO: density profile
        QS : S-wave quality factor profile
        GAMMA : incidence angle (default 0)
        POL : polarisation? (default 0)
        IDT : ?...?

    OUTPUT
        VS : complex displacement in each layer

    """
    NE = len(BETA) - 1

    A = np.zeros((2,2), dtype=complex)
    B = np.zeros((2,2), dtype=complex)
    PSH = np.zeros((2, 2, NE), dtype=complex)

    VS = np.zeros(NE+1, dtype=complex)
    ESF = np.zeros(NE+1, dtype=complex)

    UI = 1j
    UR = 1.
    ZER = 0.

    PI = np.pi
    GAMMA *= PI/180.

    POL *= PI/180.
    if np.abs(cm.cos(POL)) < 0.001:
        raise ValueError("Polarization angle too large")

    OMEGA = FREQ*2.0*PI
    if OMEGA < 0.001:
        OMEGA = 0.001

    AK = OMEGA/BETA[NE]*cm.sin(GAMMA)
    ETHS = OMEGA/BETA[NE]*cm.cos(GAMMA)
    CAMUHS = RHO[NE]*BETA[NE]**2

    for IE in range(NE):

        if QS is None or QS[IE] is None:
            CBETA = BETA[IE]
        else:
            CBETA = BETA[IE]/(1.+0.5*1j/QS[IE])

        CAMUE = RHO[IE]*CBETA**2

        ET = OMEGA/CBETA
        ET = ET*ET - AK*AK
        ET = cm.sqrt(ET)
        if np.imag(ET) < 0:
            ET = -ET
        ET = ET*H[IE]

        PSH[0][0][IE] = cm.cos(ET)
        PSH[0][1][IE] = cm.sin(ET)/ET/CAMUE*H[IE]
        PSH[1][0][IE] = -cm.sin(ET)*ET*CAMUE/H[IE]
        PSH[1][1][IE] = PSH[0][0][IE]

    A[0][0] = UR
    A[0][1] = ZER
    A[1][0] = ZER
    A[1][1] = UR

    for IE in range(NE):
        N = NE-IE-1
        for IA in [0, 1]:
            for JA in [0, 1]:
                SUM = ZER
                for K in [0, 1]:
                    SUM = SUM+A[IA][K]*PSH[K][JA][N]
                B[IA][JA] = SUM

        for IA in [0, 1]:
            for JA in [0, 1]:
                A[IA][JA] = B[IA][JA]

    if IDT == -1:
        VS[0] = 2.0/(A[0][0]+UI*A[1][0]/CAMUHS/ETHS)
    if IDT == 1:
        VS[0] = 2.0/(A[0][0]-UI*A[1][0]/CAMUHS/ETHS)

    VS[0] *= cm.cos(POL)
    DISP = VS[0]
    STRESS = ZER
    ESF[0] = ZER

    for IE in range(0, NE):
        VS[IE+1] = PSH[0][0][IE]*DISP+PSH[0][1][IE]*STRESS
        ESF[IE+1] = PSH[1][0][IE]*DISP+PSH[1][1][IE]*STRESS
        DISP = VS[IE+1]
        STRESS = ESF[IE+1]

    return VS


def haspsv(FREQ, H, ALPHA, BETA, RHO, QP, QS, IWAVE='sh', GAMMA=0., POL=90., IDT=1):
    """
    INPUT
        FREQ : calculation frequency
        H : layer thickness
        ALPHA : P-wave velocity profile
        BETA : S-wave velocity profile
        RHO: density profile
        QP : P-wave quality factor profile
        QS : S-wave quality factor profile
        IWAVE : incident wave type (1=P, 2=SV)
        GAMMA : incidence angle (default 0)
        POL : polarisation? (default 0)
        IDT : ?...?

    OUTPUT
        WS : complex VERTICAL displacement in each layer
        VS : complex HORIZONTAL displacement in each layer
    """
    NE = len(BETA) - 1

    CAUX = np.zeros(4, dtype=complex)
    DESF = np.zeros(4, dtype=complex)
    A = np.zeros((4,4), dtype=complex)
    B = np.zeros((4,4), dtype=complex)
    C = np.zeros((4,4), dtype=complex)
    E = np.zeros((4,4), dtype=complex)
    F = np.zeros((4,4), dtype=complex)
    FINV = np.zeros((4,4), dtype=complex)
    P = np.zeros((4,4, NE), dtype=complex)

    US = np.zeros(NE+1, dtype=complex)
    WS = np.zeros(NE+1, dtype=complex)

    UI = 1j
    UR = 1.
    ZER = 0.

    PI = np.pi
    GAMMA *= PI/180.

    POL *= PI/180.
    if np.abs(cm.sin(POL)) < 0.001:
        raise ValueError("Polarization angle too large")

    OMEGA = FREQ*2.0*PI
    if OMEGA < 0.001:
        OMEGA = 0.001

    OMEGA2 = OMEGA*OMEGA
    CBETHS = BETA[NE]
    CALFHS = ALPHA[NE]
    CAMUHS = CBETHS*CBETHS*RHO[NE]

    if IWAVE == 1:
        K = OMEGA*cm.sin(GAMMA)/CALFHS
    if IWAVE == 2:
        K = OMEGA*cm.sin(GAMMA)/CBETHS
    K2 = K*K

    for J in range(NE):

        if QP is None or QP[J] is None:
            CALFE = ALPHA[J]
        else:
            CALFE = ALPHA[J]/(1.+0.5*1j/QP[J])

        if QS is None or QS[J] is None:
            CBETE = BETA[J]
        else:
            CBETE = BETA[J]/(1.+0.5*1j/QS[J])

        CAMU = RHO[J]*CBETE*CBETE
        GAMA = cm.sqrt(K2-OMEGA2/CALFE/CALFE)
        NU = cm.sqrt(K2-OMEGA2/CBETE/CBETE)
        KMNU = K2+NU*NU
        FAC1 = 1.0/(OMEGA2*RHO[J])+0.0*UI
        FAC2 = FAC1*CAMU
        KAPA = GAMA*H[J]/2.0
        SEN1 = (cm.exp(KAPA)-cm.exp(-KAPA))/2.0
        SEN1 = SEN1*SEN1
        PSI = NU*H[J]/2.0
        SEN2 = (cm.exp(PSI)-cm.exp(-PSI))/2.0
        SEN2 = SEN2*SEN2
        KAPA = GAMA*H[J]
        SEN3 = (cm.exp(KAPA)-cm.exp(-KAPA))/2.0
        PSI = NU*H[J]
        SEN4 = (cm.exp(PSI)-cm.exp(-PSI))/2.0

        SUM = 2.0*K2*SEN1-KMNU*SEN2
        P[0][0][J] = 1.0+2.0*FAC2*SUM
        SUM = KMNU*SEN3/GAMA-2.0*NU*SEN4
        P[0][1][J] = K*FAC2*SUM
        SUM = K2*SEN3/GAMA-NU*SEN4
        P[0][2][J] = FAC1*SUM
        SUM = SEN1-SEN2
        P[0][3][J] = 2.0*K*FAC1*SUM
        SUM = KMNU*SEN4/NU-2.0*GAMA*SEN3
        P[1][0][J] = K*FAC2*SUM
        SUM = 2.0*K2*SEN2-KMNU*SEN1
        P[1][1][J] = 1.0+2.0*FAC2*SUM
        P[1][2][J] = -P[0][3][J]
        SUM = K2*SEN4/NU-GAMA*SEN3
        P[1][3][J] = FAC1*SUM
        SUM = 4.0*K2*GAMA*SEN3-KMNU*KMNU*SEN4/NU
        P[2][0][J] = CAMU*FAC2*SUM
        P[2][1][J] = 2.0*CAMU*CAMU*KMNU*P[0][3][J]
        P[2][2][J] = P[0][0][J]
        P[2][3][J] = -P[1][0][J]
        P[3][0][J] = -P[2][1][J]
        SUM = 4.0*K2*NU*SEN4-KMNU*KMNU*SEN3/GAMA
        P[3][1][J] = CAMU*FAC2*SUM
        P[3][2][J] = -P[0][1][J]
        P[3][3][J] = P[1][1][J]

    for I in [0, 1, 2, 3]:
        for J in [0, 1, 2, 3]:
            A[I][J] = UR
            if I != J:
                A[I][J] = ZER

    for IE in range(NE):
        N = NE-IE-1
        for IA in [0, 1, 2, 3]:
            for JA in [0, 1, 2, 3]:
                CSUM = ZER
                for KK in [0, 1, 2, 3]:
                    CSUM = CSUM+A[IA][KK]*P[KK][JA][N]
                B[IA][JA] = CSUM

        for IA in [0, 1, 2, 3]:
            for JA in [0, 1, 2, 3]:
                A[IA][JA] = B[IA][JA]

    if IWAVE == 1:
        GAMA = -UI*OMEGA*cm.cos(GAMMA)/CALFHS
        SENGAS = CBETHS*cm.sin(GAMMA)/CALFHS
        SENGAS = cm.sqrt(1.0-SENGAS*SENGAS)
        NU = -UI*OMEGA*SENGAS/CBETHS

    if IWAVE == 2:
        NU = -UI*OMEGA*cm.cos(GAMMA)/CBETHS
        SENGAS = CALFHS*cm.sin(GAMMA)/CBETHS
        SENGAS = cm.sqrt(1.0-SENGAS*SENGAS)
        GAMA = -UI*OMEGA*SENGAS/CALFHS

    KMNU = K2+NU*NU

    for IA in [0, 1, 2, 3]:
        for JA in [0, 1, 2, 3]:
            E[IA][JA] = ZER

    FAC = 2.0*CALFHS*CAMUHS*GAMA*NU*OMEGA
    FAC = CBETHS/FAC

    E[0][0] = UR*FAC
    E[1][1] = -UI*FAC
    E[2][2] = UR*FAC
    E[3][3] = -UI*FAC
    F[0][0] = 2.0*CBETHS*CAMUHS*K*GAMA*NU
    F[0][1] = -CBETHS*CAMUHS*NU*KMNU
    F[0][2] = -CBETHS*K*NU
    F[0][3] = CBETHS*GAMA*NU
    F[1][0] = -CALFHS*CAMUHS*GAMA*KMNU
    F[1][1] = 2.0*CALFHS*CAMUHS*K*GAMA*NU
    F[1][2] = CALFHS*GAMA*NU
    F[1][3] = -CALFHS*K*GAMA
    F[2][0] = 2.0*CBETHS*CAMUHS*K*GAMA*NU
    F[2][1] = CBETHS*CAMUHS*NU*KMNU
    F[2][2] = CBETHS*K*NU
    F[2][3] = CBETHS*GAMA*NU
    F[3][0] = -CALFHS*CAMUHS*GAMA*KMNU
    F[3][1] = -2.0*CALFHS*CAMUHS*K*GAMA*NU
    F[3][2] = -CALFHS*GAMA*NU
    F[3][3] = -CALFHS*K*GAMA

    for IA in [0, 1, 2, 3]:
        for JA in [0, 1, 2, 3]:
            CSUM = ZER
            for KK in [0, 1, 2, 3]:
                CSUM = CSUM+E[IA][KK]*F[KK][JA]
            FINV[IA][JA] = CSUM

    for IA in [0, 1, 2, 3]:
        for JA in [0, 1, 2, 3]:
            CSUM = ZER
            for KK in [0, 1, 2, 3]:
                CSUM = CSUM+FINV[IA][KK]*A[KK][JA]
            C[IA][JA] = CSUM

    if IWAVE == 1:
        DEN = C[3][0]*C[2][1]-C[2][0]*C[3][1]
        UZ = -C[3][1]/DEN
        WZ = C[3][0]*UI/DEN

    if IWAVE == 2:
        DEN = C[3][1]*C[2][0]-C[3][0]*C[2][1]
        UZ = -C[2][1]/DEN
        WZ = C[2][0]*UI/DEN

    US[0] = UZ*cm.sin(POL)
    WS[0] = WZ*cm.sin(POL)

    DESF[0] = US[0]
    DESF[1] = WS[0]/UI
    DESF[2] = ZER
    DESF[3] = ZER

    for IE in range(0, NE):
        for IA in [0, 1, 2, 3]:
            SUM = ZER
            for KA in [0, 1, 2, 3]:
                SUM = SUM+P[IA][KA][IE]*DESF[KA]
            CAUX[IA] = SUM

        US[IE+1] = CAUX[0]
        WS[IE+1] = CAUX[1]*UI
        DESF[0] = CAUX[0]
        DESF[1] = CAUX[1]
        DESF[2] = CAUX[2]
        DESF[3] = CAUX[3]

    return US, WS