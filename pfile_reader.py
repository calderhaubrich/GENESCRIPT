#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import sys

def ReadPfile(filepath):
    #TODO: Add check on the filepath.

    #File Converter#
    pfile = filepath.readlines()
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
    return chunked_list

#Plot Profiles Function#
#TODO: Clean up variable names.
def pfileplot(l,types,tiles,plotlabel):
    titles = l[0]
    del l[0]
    x=[]
    y=[]
    for nh in range(len(l)):
        x.append(float(l[nh][0]))
        y.append(float(l[nh][1]))
    plt.xlabel(titles[1])
    plt.ylabel(types)
    plt.plot(x,y,label=plotlabel)
    plt.legend()
    plt.title(tiles)

#Define pfile_comparison.

def main(args):
    parser = argparse.ArgumentParser(description='Reads pfile input and plots species.')
    parser.add_argument('file', type=argparse.FileType('r'), help='input some pfile')
    parser.add_argument('-l', '--list', type=int,  default=[], help='List of indices for data to plot.') #TODO: Take a list of which data to plot.
    parser.add_argument('-a', '--all',  type=bool, default=False, help='Plot everything at once.') #TODO: Just take a flag, no arguments
    #TODO: How to handle species??? Maybe not a big deal.
    args = parser.parse_args(args)

    #TODO: If no arguments passed
    #CALL pfile_comparison function.
    #else, check file exists and handle the optional flags.

    chunked_list = ReadPfile(args.file)

    if (args.all):
        fig1 = plt.figure()
        #TODO: Store these strings in some global arrays at the top, matching the indices passed in and the data indices.
        n = 'Density (n[10^20/m^3])'
        nt = 'Densities of Electrons and Ions'
        pfileplot(chunked_list[0],n,nt,'ne')
        pfileplot(chunked_list[2],n,nt,'ni')
        plt.grid() #TODO: How to handle the grid?


        fig2 = plt.figure()
        t = 'Temperature (keV)'
        tt = 'Temperatures of Electrons and Ions'
        pfileplot(chunked_list[1],t,tt,'Te')
        pfileplot(chunked_list[3],t,tt,'Ti')



        fig3 = plt.figure()
        bt = 'Density of b and nz1 Species'
        pfileplot(chunked_list[4],n,bt,'nb')
        pfileplot(chunked_list[18],n,bt,'nz1')


        fig4 = plt.figure()
        p = 'Pressure (kPa)'
        pt = 'Pressure of b Species'
        pfileplot(chunked_list[5],p,pt,'pb')


        fig5 = plt.figure()
        ppt = 'Pressures of Species'
        pfileplot(chunked_list[6],p,ppt,'ptot')


        fig6 = plt.figure()
        o = 'Angular Velocity (kRad/s)'
        ot = 'Angular Velocities'
        pfileplot(chunked_list[7],o,ot,'omeg')
        pfileplot(chunked_list[9],o,ot,'omgvb')
        pfileplot(chunked_list[11],o,ot,'omgeb')
        pfileplot(chunked_list[12],o,ot,'ommvb')



        fig7 = plt.figure()
        pfileplot(chunked_list[8],o,ot,'omegp')
        pfileplot(chunked_list[10],o,ot,'omgpp')
        pfileplot(chunked_list[13],o,ot,'ommpp')
        pfileplot(chunked_list[15],o,ot,'omepp')



        fig8 = plt.figure()
        pfileplot(chunked_list[14],o,ot,'omevb')


        fig9 = plt.figure()
        e = 'Electric Field (kV/m)'
        et = 'Electric Field'
        pfileplot(chunked_list[16],e,et,'er')


        fig10 = plt.figure()
        k = '(km/s/T)'
        pfileplot(chunked_list[17],k,k,'kpol')


        fig11 = plt.figure()
        v = 'Velocity (km/s)'
        vt = 'Toroidal Velocity'
        vp = 'Poloidal Velocity'
        pfileplot(chunked_list[19],v,vt,'vtorl')


        fig12 = plt.figure()
        pfileplot(chunked_list[20],v,vp,'vpol1')

        plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])