#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xlrd


iterpfile = '/home/calderhaubrich/Desktop/p_g_files/ITER/ITER_DT_15MA_Q10_peaked_n.txt'

def ReadPfile(filepath):
    f     = open(filepath,'r')
    pfile = f.readlines()
    del pfile[0:18]
    rawdata = []
    for i in range(len(pfile)):
        row = pfile[i].split()
        rawdata.append(row)
    profile_chunks = np.array(rawdata).T.tolist()
    return profile_chunks #profile_chunks[50] is rho_tor

def pfileplot(rho,species):
    xtitles = rho[0]
    del rho[0]
    ytitles = species[0]
    del species[0]
    x=[]
    y=[]
    for pos in range(len(rho)):
        x.append(float(rho[pos]))
    for pos2 in range(len(species)):
        y.append(float(species[pos2]))
    plt.xlabel(xtitles)
    plt.ylabel(ytitles)
    plt.plot(x,y)
    plt.grid()


pfileplot(ReadPfile(iterpfile)[50],ReadPfile(iterpfile)[15])
plt.show()
