#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.integrate as integ
import sys
import pydiag.utils.comm as com
import pydiag.utils.geom as gm
from matplotlib import ticker

from plotGrowthRates import readGrowthRates
import readPhi
import scriptUtils

phiratio        = []
ky              = []
kxgrid          = []
geometricfactor = []
def kperpSquaredCalc(phiData, geomData, common):
    phiSquared = phiData.real**2 + phiData.imag**2 #Complex square of phi.
    kymin      = common.pars.get("kymin")
    kxp        = scriptUtils.getKxArray(common)
    #Calculate kperp squared average.
    denom = 0
    numer = 0
    geom  = 0
    phiintegral = []
    for kx in range(len(kxp)):
        phiintegral.append(integ.trapz(phiSquared[kx,:]))
        denom += integ.trapz(phiSquared[kx,:])
        kxVal =  kxp[kx]
        numer += integ.trapz((geomData.gxx[:]*kxVal**2 + 2*geomData.gxy[:]*kxVal*kymin + geomData.gyy[:]*kymin**2)*phiSquared[kx,:])
        geom  += integ.trapz(geomData.gxx[:]*kxVal**2 + 2*geomData.gxy[:]*kxVal*kymin + geomData.gyy[:]*kymin**2)
    kxp, phiintegral = zip(*sorted(zip(kxp, phiintegral)))
    kperpSqAv = numer/denom
    geometricfactor.append(geom)
    phiratio.append(phiintegral/denom)
    ky.append(kymin)
    kxgrid.append(kxp)

    return kperpSqAv, denom, phiintegral

def kperpSqPlot(kperpSqAv, label):
    modes  = readGrowthRates('scan.log')
    kyRhoi = modes[0]
    gamma  = modes[2]

    plt.figure()
    plt.plot(kyRhoi, np.transpose(gamma/kperpSqAv), marker='o', color='C0')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('$\\frac{\\gamma}{\\langle k_{\\perp}^2\\rangle}$', fontsize=18)
    #plt.title('Saturation Amplitudes')
    plt.grid()

    plt.figure()
    plt.plot(kyRhoi, np.sqrt(kperpSqAv), marker='o', color='C0')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('$\\sqrt{\\langle k_{\\perp}^2\\rangle}$', fontsize=18)
    plt.grid()
    #If labels are given then plot a legend.
    if label:
        plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(kyRhoi, geometricfactor, marker='o', color='C0')
    plt.xlabel('$k_y\\rho_i$', fontsize=16)
    plt.ylabel('$\\int\\sum_{k_x}(g^{xx}k_{x}^{2} + 2g^{xy}k_{x}k_{y} + g^{yy}k_{y}^{2})dz$', fontsize=18)
    plt.grid()
    plt.tight_layout()

def main(args):
    parser = argparse.ArgumentParser("Script for calculating quasilinear saturation amplitudes.")
    parser.add_argument("directories",       nargs='+', default=[], help="Scan directories. Passing more than one shares the plot.")
    parser.add_argument("-n", "--datanames", nargs='+', default=[], help="Array of data names.")
    scriptUtils.addSavefileFlag(parser)
    fv   = parser.add_argument("-v", "--fieldVal", help="0=phi, 1=A_parallel, 2=B_parallel.", type=int, default=0)
    args = parser.parse_args(args)

    if (args.fieldVal < 0 or args.fieldVal > 2):
        parser.error('Unsupported field type: ' + str(args.fieldVal) + '. Please try one of the following: ' +
            fv.help + '.')

    combined = False
    if (len(args.directories) > 1):
        combined = True

    #If labels aren't passed in, just use blanks and no legend will be plotted.
    labels = args.datanames if len(args.datanames) > 0 else ['' for x in range(len(args.directories))]

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
        phint      = []
        denom      = []
        geomData   = gm.Geometry(common) #Geometry files are the same for all runs so just get it now.

        #Load all phi and kperp data.
        for i, extension in enumerate(extensions):
            common       = com.CommonData(extension, -1, -2)
            phiData      = readPhi.readField3D(common, extension, -1, args.fieldVal)[:, 0, :]
            kperpSqAv[i] = kperpSquaredCalc(phiData, geomData, common)[0]

        #Plot all the data.
        kygrid = np.vstack([ky]*192)
        print(np.shape(geometricfactor))
        print(kperpSqPlot(kperpSqAv, j))
        plt.figure()
        for each in range(10):
            plt.plot(kxgrid[each],phiratio[each],label=str(ky[each]))
        plt.yscale('log')
        plt.legend()
        plt.xlabel('kx')
        plt.ylabel('$\\int|\\hat{\\Phi}_{k_{x}k_{y}}(z)|^2dz$')
        plt.show()

        #contour map of phi/sumphi
        fig,ax=plt.subplots(1,1)
        #for gridspace in range(20):
        #    cp = ax.contourf(kxgrid[:][gridspace],ky, phiratio,1000,cmap = 'plasma',extend='both')
        #cp = ax.contourf(kxgrid,np.transpose(kygrid),phiratio,2000,locator=ticker.LogLocator(), cmap = 'plasma',extend='both')
        cp = ax.contourf(kxgrid,np.transpose(kygrid),phiratio,2000, cmap = 'plasma',extend='both')
        #ax.set_xlim(0,1)
        fig.colorbar(cp) # Add a colorbar to a plot
        ax.set_title('$\\int|\\hat{\\Phi}_{k_{x}k_{y}}(z)|^2dz$ \ $\\sum \\int|\\hat{\\Phi}_{k_{x}k_{y}}(z)|^2dz$')
        ax.set_xlabel('kx')
        ax.set_ylabel('ky')
        plt.show()

        #Switch back to outside scandir for saving plot.
        os.chdir(origDir)

        if (not combined):
            if (args.savefile):
                plt.savefig(args.savefile + scanDir.replace('scanfiles','') + '.pdf')
            plt.show()

    #For combining all data into one plot.
    if (combined):
        plt.grid()
        plt.legend([args.datanames[0]])
        plt.tight_layout()
        if (args.savefile):
            plt.savefig(args.savefile + 'All.pdf')
        plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
