#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys

import scriptUtils

def readGrowthRates(fileName):
    data = []
    f = open(fileName, 'r')

    for line in f.readlines():
        if (line[0] != '#'): #Ignore comment lines.
            if (len(line.split()) > 0):                 #Ignore empty lines. Not sure why there are x-values with no corresponding y sometimes.
                line = line.split()                     #Remove index from front as well.
                line[:] = [x for x in line if x != "|"] #Remove all dividers from file.
                data.append(line[1:])
    f.close()

    numKyModes = np.shape(data)[0]
    numEigVals = np.shape(data)[1]//2 #Floor by 2 to keep an integer and remove ky data.

    kyRhoi = np.zeros(numKyModes)
    gamma  = np.zeros((numEigVals, numKyModes))
    omega  = np.zeros((numEigVals, numKyModes))

    #Set up data arrays.
    for i in range(numKyModes):
        kyRhoi[i] = data[i][0]
        for j in range(numEigVals):
            gamma[j][i] = data[i][j*2 + 1] #+1 because 0 is ky modes.
            omega[j][i] = data[i][j*2 + 2]

    return [kyRhoi, omega, gamma, numEigVals]

def plotGrowthRates(title, dataPath, datafiles, datanames, norms, savefile):
    fig,axs = plt.subplots(2,1)
    kyLabel    = 'k$_y$'
    gammaLabel = '$\\gamma$'
    omegaLabel = '$\\omega$'

    if norms:
        kyLabel    += '$\\rho_i$'
        gammaLabel += ' / ($c_s$/R)'
        omegaLabel += ' / ($c_s$/R)'

    for i, datafile in enumerate(datafiles):
        [kyRhoi, omega, gamma, numEigVals] = readGrowthRates(dataPath + '/' + datafile)

        datanames[i] = datanames[i].split(',')
        for j in range(numEigVals): #Plot all eigenvalues for this data set.
            axs[0].plot(kyRhoi, gamma[j,:], marker='*', label=datanames[i][j])
            axs[1].plot(kyRhoi, omega[j,:], marker='*', label=datanames[i][j])

    plt.suptitle(title)
    axs[0].grid()
    axs[1].grid()
    axs[0].set_xlabel(kyLabel,    fontsize=14)
    axs[0].set_ylabel(gammaLabel, fontsize=14)
    axs[1].set_xlabel(kyLabel,    fontsize=14)
    axs[1].set_ylabel(omegaLabel, fontsize=14)
    axs[0].legend()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if (savefile):
        plt.savefig(dataPath + '/' + savefile + '.pdf')
    plt.show()

def main(args):
    parser = argparse.ArgumentParser("Script for plotting/saving GENE ASCII growth rates and frequencies.")
    parser.add_argument("dataPath",          help="Path to data files. Without trailing '/'.")
    parser.add_argument("-t", "--title",     default="", help="Plot title.")
    parser.add_argument("-d", "--datafiles", nargs='+', default=[], type=str, help="Array of data file names.")
    parser.add_argument("-n", "--datanames", nargs='+', default=[], type=str, help="Array of data names. If a dataset has multiple modes then separate those names with commas BUT NO spaces.")
    parser.add_argument("-z", "--norms",     action='store_true',   help="Include normalization data in labels.")
    scriptUtils.addSavefileFlag(parser)
    args = parser.parse_args(args)

    plotGrowthRates(args.title, args.dataPath, args.datafiles, args.datanames, args.norms, args.savefile)

if __name__ == "__main__":
    main(sys.argv[1:])