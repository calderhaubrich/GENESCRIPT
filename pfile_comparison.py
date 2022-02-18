#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob
import argparse
import sys



#File Converter#
unfiles = glob.glob('p163241.*')
files = sorted(unfiles,key=len)
readpfiles = []
pfilenames = list()

for file in files:
    pfilenames.append(file)
    ff = open(file,'r')
    pfile = ff.readlines()
    pfilecon = []
    chunked_list = list()
    chunk_size = 257
    for i in range(len(pfile)):
        row = pfile[i].split()
        pfilecon.append(row)
    NZA = pfilecon[4626:4630] #Species Information
    del pfilecon[4626:4630]
    for j in range(0,len(pfilecon),chunk_size):
        if len(pfilecon[j][0]) == 3:
            chunked_list.append(pfilecon[j:j+chunk_size])
        if len(pfilecon[j][0]) > 3:
            chunked_list.append(pfilecon[j:j+4])
    readpfiles.append(chunked_list)

#Pfile Comparison Plot Function
def pfileplot(l,types,name):
    plotlabel = pfilenames[name]
    del l[0]
    x=[]
    y=[]
    for nh in range(len(l)):
        x.append(float(l[nh][0]))
        y.append(float(l[nh][1]))
    plt.xlabel('psinorm')
    plt.ylabel(types)
    plt.plot(x,y,label=plotlabel)
    plt.legend()

for edens in range(len(readpfiles)):
    pfileplot(readpfiles[edens][0],'ne (10^20/m^3)',edens)
plt.show()

for etemp in range(len(readpfiles)):
    pfileplot(readpfiles[etemp][1],'te (keV)',etemp)
plt.show()

for idens in range(len(readpfiles)):
    pfileplot(readpfiles[idens][2],'ni (10^20/m^3)',idens)
plt.show()

for itemp in range(len(readpfiles)):
    pfileplot(readpfiles[itemp][3],'ti (keV)',itemp)
plt.show()

for bdens in range(len(readpfiles)):
    pfileplot(readpfiles[bdens][4],'nb (10^20/m^3)',bdens)
plt.show()

for bpres in range(len(readpfiles)):
    pfileplot(readpfiles[bpres][5],'pb (kPa)',bpres)
plt.show()

for ptot in range(len(readpfiles)):
    pfileplot(readpfiles[ptot][6],'ptot (kPa)',ptot)
plt.show()

for omeg in range(len(readpfiles)):
    pfileplot(readpfiles[omeg][7],'omeg (kRad/s)',omeg)
plt.show()

for omegp in range(len(readpfiles)):
    pfileplot(readpfiles[omegp][8],'omegp (kRad/s)',omegp)
plt.show()

for omgvb in range(len(readpfiles)):
    pfileplot(readpfiles[omgvb][9],'omgvb (kRad/s)',omgvb)
plt.show()

for omgpp in range(len(readpfiles)):
    pfileplot(readpfiles[omgpp][10],'omgpp (kRad/s)',omgpp)
plt.show()

for omgeb in range(len(readpfiles)):
    pfileplot(readpfiles[omgeb][11],'omgeb (kRad/s)',omgeb)
plt.show()

for ommvb in range(len(readpfiles)):
    pfileplot(readpfiles[ommvb][12],'omgeb (kRad/s)',omgeb)
plt.show()

for ommpp in range(len(readpfiles)):
    pfileplot(readpfiles[ommpp][13],'ommpp (kRad/s)',ommpp)
plt.show()

for omevb in range(len(readpfiles)):
    pfileplot(readpfiles[omevb][14],'omevb (kRad/s)',omevb)
plt.show()

for omepp in range(len(readpfiles)):
    pfileplot(readpfiles[omepp][15],'omepp (kRad/s)',omepp)
plt.show()

for er in range(len(readpfiles)):
    pfileplot(readpfiles[er][16],'er (kV/m)',er)
plt.show()

for kpol in range(len(readpfiles)):
    pfileplot(readpfiles[kpol][17],'kpol (km/s/T)',kpol)
plt.show()

for carb in range(len(readpfiles)):
    pfileplot(readpfiles[carb][18],'nz1 (10^20/m^3)',carb)
plt.show()

for vtor in range(len(readpfiles)):
    pfileplot(readpfiles[vtor][19],'vtor1 (km/s)',vtor)
plt.show()

for vpol in range(len(readpfiles)):
    pfileplot(readpfiles[vpol][20],'vpol1 (km/s)',vpol)
plt.show()
