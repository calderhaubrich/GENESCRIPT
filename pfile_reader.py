#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import sys
import glob
from matplotlib import cm
from tabulate import tabulate

n = 'Density (n[10^20/m^3])'
t = 'Temperature (keV)'
p = 'Pressure (kPa)'
o = 'Angular Velocity (kRad/s)'
e = 'Electric Field (kV/m)'
k = '(km/s/T)'
v = 'Velocity (km/s)'


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
def pfileplot(l,types,plotlabel):
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
    plt.grid()

#Define pfile_comparison.

def main(args):
    parser = argparse.ArgumentParser(description='Reads pfile input and plots species.')
    parser.add_argument('file', type=argparse.FileType('r'), default='p163241.03500_e0514', help='input some pfile')
    parser.add_argument('-l', '--list', nargs='*', help='List of indices for data to plot.') #TODO: Take a list of which data to plot.
    parser.add_argument('-a', '--all', default=False, action='store_true', help='Plot everything at once.')
    parser.add_argument('-c', '--comparison', default=False, action='store_true', help='Plot comparison of all pfiles.')
    args = parser.parse_args(args)

    #TODO: If no arguments passed
    #CALL pfile_comparison function.
    #else, check file exists and handle the optional flags.

    chunked_list = ReadPfile(args.file)

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

    if (args.comparison):
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

        def pfileplotcomp(l,types,name):
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

        del NZA[0][0]
        del NZA [0][3]
        del NZA [0][3]
        del NZA [0][3]

        print(tabulate(NZA))

        for edens in range(len(readpfiles)):
            pfileplotcomp(readpfiles[edens][0],'ne (10^20/m^3)',edens)
        plt.grid()
        plt.show()

        for etemp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[etemp][1],'te (keV)',etemp)
        plt.grid()
        plt.show()

        for idens in range(len(readpfiles)):
            pfileplotcomp(readpfiles[idens][2],'ni (10^20/m^3)',idens)
        plt.grid()
        plt.show()

        for itemp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[itemp][3],'ti (keV)',itemp)
        plt.grid()
        plt.show()

        for bdens in range(len(readpfiles)):
            pfileplotcomp(readpfiles[bdens][4],'nb (10^20/m^3)',bdens)
        plt.grid()
        plt.show()

        for bpres in range(len(readpfiles)):
            pfileplotcomp(readpfiles[bpres][5],'pb (kPa)',bpres)
        plt.grid()
        plt.show()

        for ptot in range(len(readpfiles)):
            pfileplotcomp(readpfiles[ptot][6],'ptot (kPa)',ptot)
        plt.grid()
        plt.show()

        for omeg in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omeg][7],'omeg (kRad/s)',omeg)
        plt.grid()
        plt.show()

        for omegp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omegp][8],'omegp (kRad/s)',omegp)
        plt.grid()
        plt.show()

        for omgvb in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omgvb][9],'omgvb (kRad/s)',omgvb)
        plt.grid()
        plt.show()

        for omgpp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omgpp][10],'omgpp (kRad/s)',omgpp)
        plt.grid()
        plt.show()

        for omgeb in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omgeb][11],'omgeb (kRad/s)',omgeb)
        plt.grid()
        plt.show()

        for ommvb in range(len(readpfiles)):
            pfileplotcomp(readpfiles[ommvb][12],'omgeb (kRad/s)',omgeb)
        plt.grid()
        plt.show()

        for ommpp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[ommpp][13],'ommpp (kRad/s)',ommpp)
        plt.grid()
        plt.show()

        for omevb in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omevb][14],'omevb (kRad/s)',omevb)
        plt.grid()
        plt.show()

        for omepp in range(len(readpfiles)):
            pfileplotcomp(readpfiles[omepp][15],'omepp (kRad/s)',omepp)
        plt.grid()
        plt.show()

        for er in range(len(readpfiles)):
            pfileplotcomp(readpfiles[er][16],'er (kV/m)',er)
        plt.grid()
        plt.show()

        for kpol in range(len(readpfiles)):
            pfileplotcomp(readpfiles[kpol][17],'kpol (km/s/T)',kpol)
        plt.grid()
        plt.show()

        for carb in range(len(readpfiles)):
            pfileplotcomp(readpfiles[carb][18],'nz1 (10^20/m^3)',carb)
        plt.grid()
        plt.show()

        for vtor in range(len(readpfiles)):
            pfileplotcomp(readpfiles[vtor][19],'vtor1 (km/s)',vtor)
        plt.grid()
        plt.show()

        for vpol in range(len(readpfiles)):
            pfileplotcomp(readpfiles[vpol][20],'vpol1 (km/s)',vpol)
        plt.grid()
        plt.show()

    #TODO: Store these strings in some global arrays at the top, matching the indices passed in and the data indices.

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
