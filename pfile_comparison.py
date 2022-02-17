#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob
import argparse
import sys



#File Converter#
files = glob.glob('p163241.*')
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


for dens in range(len(readpfiles)):
    pfileplot(readpfiles[dens][0],'ne',dens)
plt.show()
