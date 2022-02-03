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
print(chunked_list[0][0])


#Plot Profiles Density Function#
def pfileplotd(l,types):
    titles = l[0]
    del l[0]
    x=[]
    y=[]
    for n in range(len(l)):
        x.append(float(l[n][0]))
        y.append(float(l[n][1]))
    plt.xlabel(titles[1])
    plt.ylabel(types)
    plt.plot(x,y,label=titles[2])
    plt.legend()
    plt.title("Densities of Electrons and Ions")

#Plot Profiles Temperature Function#
def pfileplott(l,types):
    titles = l[0]
    del l[0]
    x=[]
    y=[]
    for n in range(len(l)):
        x.append(float(l[n][0]))
        y.append(float(l[n][1]))
    plt.xlabel(titles[1])
    plt.ylabel(types)
    plt.plot(x,y,label=titles[2])
    plt.legend()
    plt.title("Temperatures of Electrons and Ions")


n = 'Density'
pfileplotd(chunked_list[0],n)
pfileplotd(chunked_list[2],n)
plt.grid()
plt.show()

t = 'Temperature'
pfileplott(chunked_list[1],t)
pfileplott(chunked_list[3],t)
plt.grid()
plt.show()
