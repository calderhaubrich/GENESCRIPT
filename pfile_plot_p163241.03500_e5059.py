#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#File Converter#
file = open('p163241.03500_e5059', 'r')
pfile = file.readlines()
pfilecon = []
chunked_list = list()
chunk_size = 257
for i in range(len(pfile)):
    row = pfile[i].split()
    pfilecon.append(row)
NZA = pfilecon[4626:4630]
del pfilecon[4626:4630]
for j in range(0,len(pfilecon),chunk_size):
    if len(pfilecon[j][0]) == 3:
        chunked_list.append(pfilecon[j:j+chunk_size])
    if len(pfilecon[j][0]) > 3:
        chunked_list.append(pfilecon[j:j+4])

#Plot Profiles Function#
def pfileplot(l,types,tiles,plotlabel):
    titles = l[0]
    del l[0]
    x=[]
    y=[]
    for n in range(len(l)):
        x.append(float(l[n][0]))
        y.append(float(l[n][1]))
    plt.xlabel(titles[1])
    plt.ylabel(types)
    plt.plot(x,y,label=plotlabel)
    plt.legend()
    plt.title(tiles)

n = 'Density (n[10^20/m^3])'
nt = 'Densities of Electrons and Ions'
pfileplot(chunked_list[0],n,nt,'ne')
pfileplot(chunked_list[2],n,nt,'ni')
plt.grid()
plt.show()

t = 'Temperature (keV)'
tt = 'Temperatures of Electrons and Ions'
pfileplot(chunked_list[1],t,tt,'Te')
pfileplot(chunked_list[3],t,tt,'Ti')
plt.grid()
plt.show()

bt = 'Density of b and nz1 Species'
pfileplot(chunked_list[4],n,bt,'nb')
pfileplot(chunked_list[18],n,bt,'nz1')
plt.grid()
plt.show()

p = 'Pressure (kPa)'
pt = 'Pressure of b Species'
pfileplot(chunked_list[5],p,pt,'pb')
plt.grid()
plt.show()

ppt = 'Pressures of Species'
pfileplot(chunked_list[6],p,ppt,'ptot')
plt.grid()
plt.show()

o = 'Angular Velocity (kRad/s)'
ot = 'Angular Velocities'
pfileplot(chunked_list[7],o,ot,'omeg')
pfileplot(chunked_list[9],o,ot,'omgvb')
pfileplot(chunked_list[11],o,ot,'omgeb')
pfileplot(chunked_list[12],o,ot,'ommvb')
plt.grid()
plt.show()

pfileplot(chunked_list[8],o,ot,'omegp')
pfileplot(chunked_list[10],o,ot,'omgpp')
pfileplot(chunked_list[13],o,ot,'ommpp')
pfileplot(chunked_list[15],o,ot,'omepp')
plt.grid()
plt.show()

pfileplot(chunked_list[14],o,ot,'omevb')
plt.grid()
plt.show()

e = 'Electric Field (kV/m)'
et = 'Electric Field'
pfileplot(chunked_list[16],e,et,'er')
plt.grid()
plt.show()

k = '(km/s/T)'
pfileplot(chunked_list[17],k,k,'kpol')
plt.grid()
plt.show()

v = 'Velocity (km/s)'
vt = 'Toroidal Velocity'
vp = 'Poloidal Velocity'
pfileplot(chunked_list[19],v,vt,'vtorl')
plt.grid()
plt.show()

pfileplot(chunked_list[20],v,vp,'vpol1')
plt.grid()
plt.show()
