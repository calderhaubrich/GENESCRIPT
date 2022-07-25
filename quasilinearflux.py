#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
import sys
import os, fnmatch
import scipy.integrate as integ
import pydiag.utils.comm as com
import pydiag.utils.geom as gm
from plotGrowthRates import readGrowthRates
import readflux
import scriptUtils
import readPhi
import plotSatAmps

#Column Labels
datalabels = ['$\\langle|n_{1}|^{2}\\rangle$',
            '$\\langle|u_{1}|^{2}\\rangle$',
            '$\\langle|T_{\\parallel}|^{2}\\rangle$',
            '$\\langle|T_{\\perp}|^{2}\\rangle$',
            '$\\langle\\Gamma_{ES}^{x}\\rangle$',
            '$\\langle\\Gamma_{EM}^{x}\\rangle$',
            '$\\langle Q_{ES}^{x} \\rangle$',
            '$\\langle Q_{EM}^{x} \\rangle$',
            '$\\langle\\Pi_{ES}^{x}\\rangle$',
            '$\\langle\\Pi_{EM}^{x}\\rangle$']

def main(args):
    parser = argparse.ArgumentParser(description="Script for plotting/saving GENE ASCII flux output.")
    parser.add_argument("directories",       nargs='+', default=[], help="Scan directories. Passing more than one shares the plot.")
    fv   = parser.add_argument("-v", "--fieldVal", help="0=phi, 1=A_parallel, 2=B_parallel.", type=int, default=0)
    args = parser.parse_args(args)

    for j, scanDir in enumerate(args.directories):
        os.chdir(scanDir)
        extensions = scriptUtils.getFileExtensions()
        common     = com.CommonData(extensions[0], -1, -2)
        nx0        = common.pars.get("nx0")
        nz         = common.pars.get("nz0")
        kySize     = len(extensions)
        kperpSqAv  = np.zeros(kySize)
        geomData   = gm.Geometry(common)

    modes  = readGrowthRates('scan.log')
    kyRhoi = modes[0]
    gamma  = modes[2]

    denom  = [] #TODO: This should all get simpler with nested loops and arrays.
    e_heat = []
    i_heat = []
    w_heat = []
    e_flux = []
    i_flux = []
    w_flux = []
    for i, extension in enumerate(extensions):
        common       = com.CommonData(extension, -1, -2)
        nrgpath      = 'nrg' + extension
        nrgdata      = os.getcwd() +'/'+ nrgpath
        fluxvalues = readflux.readFlux(nrgdata)
        e_heat.append(fluxvalues[0][-1])
        i_heat.append(fluxvalues[1][-1])
        w_heat.append(fluxvalues[2][-1])
        e_flux.append(fluxvalues[3][-1])
        i_flux.append(fluxvalues[4][-1])
        w_flux.append(fluxvalues[5][-1])
        phiData      = readPhi.readField3D(common, extension, -1, args.fieldVal)[:, 0, :]
        phidenom     = readPhi.readField3D(common, extension, -1, args.fieldVal)[0,0,0]
        kperpSqAv[i] = plotSatAmps.kperpSquaredCalc(phiData, geomData, common)[0]
        normilizer = phidenom.real**2 + phidenom.imag**2
        denom.append(normilizer**2)
    ampshape   = (np.transpose(gamma/kperpSqAv))

    qiheat = []
    qeheat = []
    qwheat = []
    qiflux = []
    qeflux = []
    qwflux = []
    ratio  = []
    for ql in range(len(ampshape)):
        quasi_i_heat = (ampshape[ql]/denom[ql]) * float(i_heat[ql])
        qiheat.append(quasi_i_heat)
        quasi_e_heat = (ampshape[ql]/denom[ql]) * float(e_heat[ql])
        qeheat.append(quasi_e_heat)
        quasi_w_heat = (ampshape[ql]/denom[ql]) * float(w_heat[ql])
        qwheat.append(quasi_w_heat)
        quasi_i_flux = (ampshape[ql]/denom[ql]) * float(i_flux[ql])
        qiflux.append(quasi_i_flux)
        quasi_e_flux = (ampshape[ql]/denom[ql]) * float(e_flux[ql])
        qeflux.append(quasi_e_flux)
        quasi_w_flux = (ampshape[ql]/denom[ql]) * float(w_flux[ql])
        qwflux.append(quasi_w_flux)
        ratio.append(quasi_i_flux/quasi_i_heat)
    
    print('Q.L. ion heat  = {}'.format(readflux.quasisum(qiheat)))
    print('Q.L. ele. heat = {}'.format(readflux.quasisum(qeheat)))
    print('Q.L. W heat    = {}'.format(readflux.quasisum(qwheat)))
    print('Q.L. ion flux  = {}'.format(readflux.quasisum(qiflux)))
    print('Q.L. ele. flux = {}'.format(readflux.quasisum(qeflux)))
    print('Q.L. W flux    = {}'.format(readflux.quasisum(qwflux)))

    plt.figure() #TODO: What I've been doing is just having main call functions like "getData()" and "plotData()" essentially. So these scripts are easy to read.
    plt.plot(kyRhoi, qiheat, label = 'Q$_i$')
    plt.plot(kyRhoi, qeheat, label = 'Q$_e$')
    plt.plot(kyRhoi, qwheat, label = 'Q$_W$')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('Q$_{QL}$', fontsize=18)
    #plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.plot(kyRhoi, qiflux, label = '$\\Gamma_{i}$')
    plt.plot(kyRhoi, qeflux, label = '$\\Gamma_{e}$')
    plt.plot(kyRhoi, qwflux, label = '$\\Gamma_{W}$')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('$\\Gamma_{QL}$', fontsize=18)
    #plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    plt.figure()
    plt.plot(kyRhoi, ratio, label = '$\\Gamma_{i}$ / Q$_i$')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('$\\Gamma_{i}$ / Q$_i$', fontsize=18)
    plt.grid()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
