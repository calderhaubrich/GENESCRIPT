#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import sys
import glob
from tabulate import tabulate

n = 'Density (n[10^20/m^3])'
t = 'Temperature (keV)'
p = 'Pressure (kPa)'
o = 'Angular Velocity (kRad/s)'
e = 'Electric Field (kV/m)'
k = '(km/s/T)'
v = 'Velocity (km/s)'

colors = ['darkblue','mediumblue','blue','dodgerblue', 'mediumturquoise','aquamarine','mediumspringgreen','palegreen','chartreuse','greenyellow','yellow']

#Individual pfile reader for input file#
def ReadPfile(filepath):
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
def pfileplot(heads,types,plotlabel):
    titles = heads[0]
    del heads[0]
    x=[]
    y=[]
    for pos in range(len(heads)):
        x.append(float(heads[pos][0]))
        y.append(float(heads[pos][1]))
    plt.xlabel("$\\rho_{tor}$", fontsize=16)
    plt.ylabel(types, fontsize=16)
    for e in range(len(x)):
        y[e] = y[e] * (10**20)
    plt.plot(x,y,label=plotlabel,linewidth=3.0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend()
    plt.grid()

def main(args):
    parser = argparse.ArgumentParser(description='Reads pfile input and plots species.')
    parser.add_argument('file', type=argparse.FileType('r'), default='p163241.03500_e0514', help='input some pfile')
    parser.add_argument('-l', '--list', nargs='*', help='List of indices for data to plot.')
    parser.add_argument('-a', '--all', default=False, action='store_true', help='Plot everything at once.')
    parser.add_argument('-c', '--comparison', default=False, action='store_true', help='Plot comparison of all pfiles.')
    parser.add_argument('-i','--info', default=False, action='store_true', help='gives documentation on plots')
    args = parser.parse_args(args)

    if (args.info):
        print("""
        This code plots pfile species for case 163241.03500.

        Input: -a or --all
        Plots all species
        
        Input: -c or --comparison
        Plots all species together over all time stamps

        Input: -l or --list
        Input list data: 0 1 2 3 4 etc...
        0  = ne
        1  = te
        2  = ni
        3  = ti
        4  = nb
        5  = pb
        6  = ptot
        7  = omeg
        8  = omegp
        9  = omgvb
        10 = omgpp
        11 = omgeb
        12 = ommvb
        13 = ommpp
        14 = omevb
        15 = omepp
        16 = er
        17 = kpol
        18 = nzl
        19 = vtor1
        20 = vpol1
        """)

    chunked_list = ReadPfile(args.file)

    #Plots all species with respect to psinorm#
    if (args.all):
        pfileplot(chunked_list[0],n,'ne')
        plt.show()

        pfileplot(chunked_list[1],t,'te')
        plt.show()

        pfileplot(chunked_list[2],n,'ni')
        plt.show()

        pfileplot(chunked_list[3],t,'ti')
        plt.show()

        pfileplot(chunked_list[4],n,'nb')
        plt.show()

        pfileplot(chunked_list[5],p,'pb')
        plt.show()

        pfileplot(chunked_list[6],p,'ptot')
        plt.show()

        pfileplot(chunked_list[7],o,'omeg')
        plt.show()

        pfileplot(chunked_list[8],o,'omegp')
        plt.show()

        pfileplot(chunked_list[9],o,'omgvb')
        plt.show()

        pfileplot(chunked_list[10],o,'omgpp')
        plt.show()

        pfileplot(chunked_list[11],o,'omgeb')
        plt.show()

        pfileplot(chunked_list[12],o,'ommvb')
        plt.show()

        pfileplot(chunked_list[13],o,'ommpp')
        plt.show()

        pfileplot(chunked_list[14],o,'omevb')
        plt.show()

        pfileplot(chunked_list[15],o,'omepp')
        plt.show()

        pfileplot(chunked_list[16],e,'er')
        plt.show()

        pfileplot(chunked_list[17],k,'kpol')
        plt.show()

        pfileplot(chunked_list[18],n,'nz1')
        plt.show()

        pfileplot(chunked_list[19],v,'vtorl')
        plt.show()

        pfileplot(chunked_list[20],v,'vpol1')
        plt.show()
        
        pfileplot(chunked_list[0],n,'ne')
        pfileplot(chunked_list[2],n,'ni')
        plt.grid()
        plt.legend()
        plt.show()

        pfileplot(chunked_list[1],t,'Te')
        pfileplot(chunked_list[3],t,'Ti')
        plt.grid()
        plt.legend()
        plt.show()

        pfileplot(chunked_list[0],n,'ne')
        pfileplot(chunked_list[2],n,'ni')
        pfileplot(chunked_list[18],n,'nz1')
        plt.grid()
        plt.legend()
        plt.show()

    #Looks for all pfiles in file and plots time comparison by species#
    if (args.comparison):
        unfiles = glob.glob('/home/calderhaubrich/Desktop/p_g_files/DIII-D/163241.03500/p163241.*')
        files = sorted(unfiles)
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

        def pfileplotcomp(l,types,name,col):
            plotlabel = pfilenames[name]
            del l[0]
            x=[]
            y=[]
            for nh in range(len(l)):
                x.append(float(l[nh][0]))
                y.append(float(l[nh][1]))
            plt.xlabel('psinorm',fontsize=16)
            plt.ylabel(types,fontsize=16)
            plt.plot(x,y,color=col,label=plotlabel)
            #plt.legend()

        for edens in range(len(readpfiles)):
            c = colors[edens]
            pfileplotcomp(readpfiles[edens][0],'n$_e$ (10$^{20}$/m$^{3}$)',edens,c)
        plt.grid()
        plt.show()

        for etemp in range(len(readpfiles)):
            c = colors[etemp]
            pfileplotcomp(readpfiles[etemp][1],'T$_e$ (keV)',etemp,c)
        plt.grid()
        plt.show()

        for idens in range(len(readpfiles)):
            c = colors[idens]
            pfileplotcomp(readpfiles[idens][2],'n$_i$ (10$^{20}$/m$^{3}$)',idens,c)
        plt.grid()
        plt.show()

        for itemp in range(len(readpfiles)):
            c = colors[itemp]
            pfileplotcomp(readpfiles[itemp][3],'T$_i$ (keV)',itemp,c)
        plt.grid()
        plt.show()

        for bdens in range(len(readpfiles)):
            c = colors[bdens]
            pfileplotcomp(readpfiles[bdens][4],'nb (10^20/m^3)',bdens,c)
        plt.grid()
        plt.show()

        for bpres in range(len(readpfiles)):
            c = colors[bpres]
            pfileplotcomp(readpfiles[bpres][5],'pb (kPa)',bpres,c)
        plt.grid()
        plt.show()

        for ptot in range(len(readpfiles)):
            c = colors[ptot]
            pfileplotcomp(readpfiles[ptot][6],'ptot (kPa)',ptot,c)
        plt.grid()
        plt.show()

        for omeg in range(len(readpfiles)):
            c = colors[omeg]
            pfileplotcomp(readpfiles[omeg][7],'omeg (kRad/s)',omeg,c)
        plt.grid()
        plt.show()

        for omegp in range(len(readpfiles)):
            c = colors[omegp]
            pfileplotcomp(readpfiles[omegp][8],'omegp (kRad/s)',omegp,c)
        plt.grid()
        plt.show()

        for omgvb in range(len(readpfiles)):
            c = colors[omgvb]
            pfileplotcomp(readpfiles[omgvb][9],'omgvb (kRad/s)',omgvb,c)
        plt.grid()
        plt.show()

        for omgpp in range(len(readpfiles)):
            c = colors[omgpp]
            pfileplotcomp(readpfiles[omgpp][10],'omgpp (kRad/s)',omgpp,c)
        plt.grid()
        plt.show()

        for omgeb in range(len(readpfiles)):
            c = colors[omgeb]
            pfileplotcomp(readpfiles[omgeb][11],'omgeb (kRad/s)',omgeb,c)
        plt.grid()
        plt.show()

        for ommvb in range(len(readpfiles)):
            c = colors[ommvb]
            pfileplotcomp(readpfiles[ommvb][12],'omgeb (kRad/s)',omgeb,c)
        plt.grid()
        plt.show()

        for ommpp in range(len(readpfiles)):
            c = colors[ommpp]
            pfileplotcomp(readpfiles[ommpp][13],'ommpp (kRad/s)',ommpp,c)
        plt.grid()
        plt.show()

        for omevb in range(len(readpfiles)):
            c = colors[omevb]
            pfileplotcomp(readpfiles[omevb][14],'omevb (kRad/s)',omevb,c)
        plt.grid()
        plt.show()

        for omepp in range(len(readpfiles)):
            c = colors[omepp]
            pfileplotcomp(readpfiles[omepp][15],'omepp (kRad/s)',omepp,c)
        plt.grid()
        plt.show()

        for er in range(len(readpfiles)):
            c = colors[er]
            pfileplotcomp(readpfiles[er][16],'er (kV/m)',er,c)
        plt.grid()
        plt.show()

        for kpol in range(len(readpfiles)):
            c = colors[kpol]
            pfileplotcomp(readpfiles[kpol][17],'kpol (km/s/T)',kpol,c)
        plt.grid()
        plt.show()

        for carb in range(len(readpfiles)):
            c = colors[carb]
            pfileplotcomp(readpfiles[carb][18],'n$_C$ (10$^{20}$/m$^{3}$)',carb,c)
        plt.grid()
        plt.show()

        for vtor in range(len(readpfiles)):
            c = colors[vtor]
            pfileplotcomp(readpfiles[vtor][19],'vtor1 (km/s)',vtor,c)
        plt.grid()
        plt.show()

        for vpol in range(len(readpfiles)):
            c = colors[vpol]
            pfileplotcomp(readpfiles[vpol][20],'vpol1 (km/s)',vpol,c)
        plt.grid()
        plt.show()

    #Input number for individual species plots with respect to psinorm#
    if (args.list):
        lists = args.list
        for element in lists:
            if float(element) == 0:
                listlabel = n
                pfileplot(chunked_list[int(element)],listlabel,'ne')
                plt.show()
            if float(element) == 1:
                listlabel = t
                pfileplot(chunked_list[int(element)],listlabel,'te')
                plt.show()
            if float(element) == 2:
                listlabel = n
                pfileplot(chunked_list[int(element)],listlabel,'ni')
                plt.show()
            if float(element) == 3:
                listlabel = t
                pfileplot(chunked_list[int(element)],listlabel,'ti')
                plt.show()
            if float(element) == 4:
                listlabel = n
                pfileplot(chunked_list[int(element)],listlabel,'nb')
                plt.show()
            if float(element) == 5:
                listlabel = p
                pfileplot(chunked_list[int(element)],listlabel,'pb')
                plt.show()
            if float(element) == 6:
                listlabel = p
                pfileplot(chunked_list[int(element)],listlabel,'ptot')
                plt.show()
            if float(element) == 7:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omeg')
                plt.show()
            if float(element) == 8:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omegp')
                plt.show()
            if float(element) == 9:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omgvb')
                plt.show()
            if float(element) == 10:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omgpp')
                plt.show()
            if float(element) == 11:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omgeb')
                plt.show()
            if float(element) == 12:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'ommvb')
                plt.show()
            if float(element) == 13:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'ommpp')
                plt.show()
            if float(element) == 14:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omevb')
                plt.show()
            if float(element) == 15:
                listlabel = o
                pfileplot(chunked_list[int(element)],listlabel,'omepp')
                plt.show()
            if float(element) == 16:
                listlabel = e
                pfileplot(chunked_list[int(element)],listlabel,'er')
                plt.show()
            if float(element) == 17:
                listlabel = k
                pfileplot(chunked_list[int(element)],listlabel,'kpol')
                plt.show()
            if float(element) == 18:
                listlabel = n
                pfileplot(chunked_list[int(element)],listlabel,'nzl')
                plt.tight_layout()
                plt.show()
            if float(element) == 19:
                listlabel = v
                pfileplot(chunked_list[int(element)],listlabel,'vtor1')
                plt.show()
            if float(element) == 20:
                listlabel = v
                pfileplot(chunked_list[int(element)],listlabel,'vpol1')
                plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
