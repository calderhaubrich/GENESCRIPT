#!/usr/bin/env python
import os
import math
import numpy as np

def quasisum(array):
    value = 0
    for v in range(len(array)):
        if math.isnan(array[v]):
            array[v] = 0
        value += array[v]
    return value[0]

def readFlux(fileName):
    f = open(fileName, 'r')
    f = f.readlines()
    nrg = []
    for line in range(len(f)):
        row = f[line].split()
        nrg.append(row)
    groups = np.array_split(nrg, int(len(f)/4))

    t          = [groups[timestamp][0][0] for timestamp in range(len(groups))]
    ion_n      = [groups[ionn][1][0] for ionn in range(len(groups))]
    ion_u      = [groups[ionu][1][1] for ionu in range(len(groups))]
    ion_T_para = [groups[iontpara][1][2] for iontpara in range(len(groups))]
    ion_T_perp = [groups[iontperp][1][3] for iontperp in range(len(groups))]
    ion_gam_es = [groups[ionges][1][4] for ionges in range(len(groups))]
    ion_gam_em = [groups[iongem][1][5] for iongem in range(len(groups))]
    ion_q_es   = [groups[ionqes][1][6] for ionqes in range(len(groups))]
    ion_q_em   = [groups[ionqem][1][7] for ionqem in range(len(groups))]
    ion_pi_es  = [groups[ionpis][1][8] for ionpis in range(len(groups))]
    ion_pi_em  = [groups[ionpim][1][9] for ionpim in range(len(groups))]
    e_n        = [groups[en][2][0] for en in range(len(groups))]
    e_u        = [groups[eu][2][1] for eu in range(len(groups))]
    e_T_para   = [groups[etpara][2][2] for etpara in range(len(groups))]
    e_T_perp   = [groups[etperp][2][3] for etperp in range(len(groups))]
    e_gam_es   = [groups[eges][2][4] for eges in range(len(groups))]
    e_gam_em   = [groups[egem][2][5] for egem in range(len(groups))]
    e_q_es     = [groups[eqes][2][6] for eqes in range(len(groups))]
    e_q_em     = [groups[eqem][2][7] for eqem in range(len(groups))]
    e_pi_es    = [groups[epis][2][8] for epis in range(len(groups))]
    e_pi_em    = [groups[epim][2][9] for epim in range(len(groups))]
    w_n        = [groups[wn][3][0] for wn in range(len(groups))]
    w_u        = [groups[wu][3][1] for wu in range(len(groups))]
    w_T_para   = [groups[wtpara][3][2] for wtpara in range(len(groups))]
    w_T_perp   = [groups[wtperp][3][3] for wtperp in range(len(groups))]
    w_gam_es   = [groups[wges][3][4] for wges in range(len(groups))]
    w_gam_em   = [groups[wgem][3][5] for wgem in range(len(groups))]
    w_q_es     = [groups[wqes][3][6] for wqes in range(len(groups))]
    w_q_em     = [groups[wqem][3][7] for wqem in range(len(groups))]
    w_pi_es    = [groups[wpis][3][8] for wpis in range(len(groups))]
    w_pi_em    = [groups[wpim][3][9] for wpim in range(len(groups))]
    return e_q_es, ion_q_es, w_q_es, e_gam_es, ion_gam_es,w_gam_es
