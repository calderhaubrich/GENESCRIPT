#!/usr/bin/env python

import argparse
#from re import I
#from this import d
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
import scriptUtils
import readPhi
import plotSatAmps

#Column Labels
n      = '$\\langle|n_{1}|^{2}\\rangle$' #TODO: Use an array.
u      = '$\\langle|u_{1}|^{2}\\rangle$'
T_para = '$\\langle|T_{\\parallel}|^{2}\\rangle$'
T_perp = '$\\langle|T_{\\perp}|^{2}\\rangle$'
Gam_es = '$\\langle\\Gamma_{ES}^{x}\\rangle$'
Gam_em = '$\\langle\\Gamma_{EM}^{x}\\rangle$'
q_es   = '$\\langle Q_{ES}^{x} \\rangle$'
q_em   = '$\\langle Q_{EM}^{x} \\rangle$'
pi_es  = '$\\langle\\Pi_{ES}^{x}\\rangle$'
pi_em  = '$\\langle\\Pi_{EM}^{x}\\rangle$'

fluxLabels = [n,u,T_para,T_perp,Gam_es,Gam_em,q_es,q_em,pi_es,pi_em]

t      = [] #time
i_data = [] #ion information
e_data = [] #electron information
w_data = [] #tungsten information

#ion data groups #Use arrays of arrays, like [species][value]
ion_n      = [] 
ion_u      = []
ion_T_para = []
ion_T_perp = []
ion_gam_es = []
ion_gam_em = []
ion_q_es   = []
ion_q_em   = []
ion_pi_es  = []
ion_pi_em  = []

#electron data groups
e_n        = []
e_u        = []
e_T_para   = []
e_T_perp   = []
e_gam_es   = []
e_gam_em   = []
e_q_es     = []
e_q_em     = []
e_pi_es    = []
e_pi_em    = []

#tungsten data groups
w_n        = []
w_u        = []
w_T_para   = []
w_T_perp   = []
w_gam_es   = []
w_gam_em   = []
w_q_es     = []
w_q_em     = []
w_pi_es    = []
w_pi_em    = []

def find(name, path):
    result = []
    for root, dirs, fies in os.walk(path):
        if name in fies:
            return os.path.join(root,name)

def readFlux(fileName): #TODO: This should be in a separate script for reading/plotting nrg file data. I have a file but cant recall if it reads the nrg or gui output.
    f = open(fileName, 'r')
    f = f.readlines()
    nrg = []
    for line in range(len(f)):
        row = f[line].split()
        nrg.append(row)
    dividor = int(len(f)/4)
    groups = np.array_split(nrg, dividor)
    for timestamp in range(len(groups)):    #TODO: This stuff should agnostically loop over n species. Probably can be set up depending on the file structure.
        t.append(float(groups[timestamp][0][0]))
        i_data.append(groups[timestamp][1])
        e_data.append(groups[timestamp][2])
        w_data.append(groups[timestamp][3])
    for datapoint in range(len(t)):
        ion_n.append(float(i_data[datapoint][0]))
        ion_u.append(float(i_data[datapoint][1]))
        ion_T_para.append(float(i_data[datapoint][2]))
        ion_T_perp.append(float(i_data[datapoint][3]))
        ion_gam_es.append(float(i_data[datapoint][4]))
        ion_gam_em.append(float(i_data[datapoint][5]))
        ion_q_es.append(float(i_data[datapoint][6]))
        ion_q_em.append(float(i_data[datapoint][7]))
        ion_pi_es.append(float(i_data[datapoint][8]))
        ion_pi_em.append(float(i_data[datapoint][9]))
        e_n.append(float(e_data[datapoint][0]))
        e_u.append(float(e_data[datapoint][1]))
        e_T_para.append(float(e_data[datapoint][2]))
        e_T_perp.append(float(e_data[datapoint][3]))
        e_gam_es.append(float(e_data[datapoint][4]))
        e_gam_em.append(float(e_data[datapoint][5]))
        e_q_es.append(float(e_data[datapoint][6]))
        e_q_em.append(float(e_data[datapoint][7]))
        e_pi_es.append(float(e_data[datapoint][8]))
        e_pi_em.append(float(e_data[datapoint][9]))
        w_n.append(float(w_data[datapoint][0]))
        w_u.append(float(w_data[datapoint][1]))
        w_T_para.append(float(w_data[datapoint][2]))
        w_T_perp.append(float(w_data[datapoint][3]))
        w_gam_es.append(float(w_data[datapoint][4]))
        w_gam_em.append(float(w_data[datapoint][5]))
        w_q_es.append(float(w_data[datapoint][6]))
        w_q_em.append(float(w_data[datapoint][7]))
        w_pi_es.append(float(w_data[datapoint][8]))
        w_pi_em.append(float(w_data[datapoint][9]))
    return e_q_es, ion_q_es, w_q_es, e_gam_es, ion_gam_es,w_gam_es

def main(args):
    parser = argparse.ArgumentParser(description="Script for plotting/saving GENE ASCII flux output.")
    parser.add_argument("directories",       nargs='+', default=[], help="Scan directories. Passing more than one shares the plot.")
    parser.add_argument('-q','--quasi', default=False, action='store_true', help='Does Quasi-Linear Analysis of flux data.')
    fv   = parser.add_argument("-v", "--fieldVal", help="0=phi, 1=A_parallel, 2=B_parallel.", type=int, default=0)
    args = parser.parse_args(args)

    #Quasi-Linear Study  #TODO: We should make this the main function of this script and always have it happen. Once flux plotting is separated out.
    if (args.quasi):     #TODO: Then we can make the input arguments more elaborate here for different quasilinear interests. Neeraj liked to plot EM/ES ratios for instance.
        path = args.directories
        for j, scanDir in enumerate(args.directories):
            #Make it so this script can be called from above a scan dir full of files.
            origDir = os.getcwd()
            os.chdir(scanDir)
            extensions = scriptUtils.getFileExtensions()
            #Set up necessary starting data w/ first param file.
            common     = com.CommonData(extensions[0], -1, -2)
            nx0        = common.pars.get("nx0")
            nz         = common.pars.get("nz0")
            kySize     = len(extensions)
            kperpSqAv  = np.zeros(kySize)
            geomData   = gm.Geometry(common) #Geometry files are the same for all runs so just get it now.

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
            fluxvalues   = readFlux(nrgdata)
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
            quasi_i_heat = (ampshape[ql]/denom[ql]) * i_heat[ql]
            qiheat.append(quasi_i_heat)
            quasi_e_heat = (ampshape[ql]/denom[ql]) * e_heat[ql]
            qeheat.append(quasi_e_heat)
            quasi_w_heat = (ampshape[ql]/denom[ql]) * w_heat[ql]
            qwheat.append(quasi_w_heat)
            quasi_i_flux = (ampshape[ql]/denom[ql]) * i_flux[ql]
            qiflux.append(quasi_i_flux)
            quasi_e_flux = (ampshape[ql]/denom[ql]) * e_flux[ql]
            qeflux.append(quasi_e_flux)
            quasi_w_flux = (ampshape[ql]/denom[ql]) * w_flux[ql]
            qwflux.append(quasi_w_flux)
            ratio.append(quasi_i_flux/quasi_i_heat)
        
        print('Q.L. ion heat  = {}'.format(np.sum(qiheat)))
        print('Q.L. ele. heat = {}'.format(np.sum(qeheat)))
        print('Q.L. W heat    = {}'.format(np.sum(qwheat)))
        print('Q.L. ion flux  = {}'.format(np.sum(qiflux)))
        print('Q.L. ele. flux = {}'.format(np.sum(qeflux)))
        print('Q.L. W flux    = {}'.format(np.sum(qwflux)))

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
