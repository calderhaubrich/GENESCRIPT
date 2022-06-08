#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cheasefiles
import millergeometryfunction
import rho_tor_to_miller_r as miller
import W_ionization_function as ionize
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from scipy.interpolate import InterpolatedUnivariateSpline, make_interp_spline, BSpline

plt.rcParams.update({'font.size':14})

proffpath = '/home/calderhaubrich/GENE/DIIID_163241_profiles/p163241.03500_e5059'

def main(args):
    parser = argparse.ArgumentParser(description='Calculates gradients at given rhotor position.')
    parser.add_argument('file', help='input some gfile')
    parser.add_argument('location', help='rhotor location')
    parser.add_argument('-r','--rho', default=False, action='store_true', help='Finds rho gradients')
    parser.add_argument('-s','--rho_star',default=False, action='store_true', help='Applies scaling factor to profiles')
    parser.add_argument('-m','--miller',default=False, action='store_true', help='Plots with miller r')
    parser.add_argument('-a','--all', default=False, action='store_true', help='Plot species with respect to rhotor')
    parser.add_argument('-i','--info', default=False, action='store_true', help='Provides info on code')
    args = parser.parse_args(args)

    efitfpath = args.file
    rho_tr = float(args.location)

    setparam = {'nrhomesh':'rhotor'}

    efitprof = cheasefiles.read_eqdsk(eqdskfpath=efitfpath,setParam=setparam)
    specprof = cheasefiles.read_profiles(profilesfpath=proffpath,setParam=setparam,eqdsk=efitfpath)


    rhotor= np.array(specprof['rhotor'])

    Te= np.array(specprof['Te'])
    Ti= np.array(specprof['Ti'])
    ne= np.array(specprof['ne'])
    ni= np.array(specprof['ni'])
    nb= np.array(specprof['nb'])
    nc= np.array(specprof['nz'])
    Pb= np.array(specprof['Pb'])
    Vtor= np.array(specprof['Vtor'])
    Vpol= np.array(specprof['Vpol'])
    Ptot= np.array(specprof['pressure'])
    omeg = np.array(specprof['omeg'])
    q= np.array(efitprof['q'])
    Tb= (Pb/nb)*(1/(1.6*10e-19))

    geomdata = millergeometryfunction.finder(efitfpath,rho_tr)
    R0 = geomdata[0]
    a  = geomdata[1]
    drho_tor_dr = geomdata[2]
    Bref = geomdata[8]

    geomiller = miller.findmiller(efitfpath)
    Bref_tr   = geomiller[2]
    Lref_tr   = geomiller[3]

    interpol_order = 3

    ne_spl = InterpolatedUnivariateSpline(rhotor, ne, k=interpol_order)
    ni_spl = InterpolatedUnivariateSpline(rhotor, ni, k=interpol_order)
    nb_spl = InterpolatedUnivariateSpline(rhotor, nb, k=interpol_order)
    te_spl = InterpolatedUnivariateSpline(rhotor, Te, k=interpol_order)
    ti_spl = InterpolatedUnivariateSpline(rhotor, Ti, k=interpol_order)
    tb_spl = InterpolatedUnivariateSpline(rhotor, Tb, k=interpol_order)
    pb_spl = InterpolatedUnivariateSpline(rhotor, Pb, k=interpol_order)
    vtor_spl = InterpolatedUnivariateSpline(rhotor, Vtor, k=interpol_order)
    Ptot_spl = InterpolatedUnivariateSpline(rhotor, Ptot, k=interpol_order)
    omeg_spl = InterpolatedUnivariateSpline(rhotor, omeg, k=interpol_order)

    dnedrho_spl = ne_spl.derivative()
    dnidrho_spl = ni_spl.derivative()
    dnbdrho_spl = nb_spl.derivative()
    dtedrho_spl = te_spl.derivative()
    dtidrho_spl = ti_spl.derivative()
    dtbdrho_spl = tb_spl.derivative()
    dpbdrho_spl = pb_spl.derivative()
    dvtordrho_spl = vtor_spl.derivative()   

    #dvpolrho = dvpoldrho_spl(rhotor)


    # Normalized paramters at target toroidal radial location
    Z0= 0.0031727/R0
    mD= 3.343583719e-27         # Mass of deuterium in kg
    m_ref= mD                   # Deuteriam mass as reference mass
    mD_tr= mD/m_ref             # Normalized deuterium mass
    me= 9.10938356e-31          # Mass of electron in SI Units (kg)
    me_tr= me/m_ref             # Normalized electron mass
    #n_ref= ne_spl(rhotor)      # reference temparture as function of rhotor
    t_ref= ti_spl(rhotor) 
    nref_tr= ne_spl(rho_tr)     # Electron density as refernce density at target location
    tref_tr= ti_spl(rho_tr)     # main ion tempeartrure as reference temperature at target location

    vtor_tr= vtor_spl(rho_tr)   # toroidal rotation velocity m/s at reference location
    omeg_tr= omeg_spl(rho_tr)   # toroidal rotation frequency in rad/sec at reference location
    ne_tr= ne_spl(rho_tr)/nref_tr
    ni_tr= ni_spl(rho_tr)/nref_tr
    nb_tr= nb_spl(rho_tr)/nref_tr

    te_tr= te_spl(rho_tr)/tref_tr
    ti_tr= ti_spl(rho_tr)/tref_tr
    tb_tr= tb_spl(rho_tr)/tref_tr


    #Calculate Gradients with respect to major radius R0
    omne_tr= dnedrho_spl(rho_tr)*[-R0/ne_spl(rho_tr)]*drho_tor_dr #omne_tr = -(Lref/n)(dn/drho)
    omni_tr= dnidrho_spl(rho_tr)*[-R0/ni_spl(rho_tr)]*drho_tor_dr
    omnb_tr= dnbdrho_spl(rho_tr)*[-R0/nb_spl(rho_tr)]*drho_tor_dr
    omte_tr= dtedrho_spl(rho_tr)*[-R0/te_spl(rho_tr)]*drho_tor_dr
    omti_tr= dtidrho_spl(rho_tr)*[-R0/ti_spl(rho_tr)]*drho_tor_dr
    omtb_tr= dtbdrho_spl(rho_tr)*[-R0/tb_spl(rho_tr)]*drho_tor_dr

    # To calculate ion density when no Carbon impurity present
    #nc_tr= 0
    #omnc_tr=0
    #ni_tr=(ne_tr-6*nc_tr);
    #omni_tr= (ne_tr*omne_tr- 6*nc_tr*omnc_tr)/ni_tr;

    # Carboon density profile from the quasineutrality condition
    nc_tr=(ne_tr - ni_tr)/6
    omnc_tr =(ne_tr*omne_tr -ni_tr*omni_tr)/(6*nc_tr)

    rho_star = (m_ref*np.sqrt(tref_tr/m_ref))/(Bref*R0) #Gyroradius-to-machine-size ratio at reference location
    rho_star_tr = (m_ref*np.sqrt(t_ref/m_ref))/(Bref_tr*Lref_tr)

    #Code Information and Short Cuts#
    if (args.info):
        print("""
        This code is used to calculate gradients with respect to major radius, plot species with respect to rhotor,
        and plots miller_r relationship.

        Input: -r or --rho
        Returns gradient values with respect to major radius R0.

        Input: -a or --all
        Returns plots of species with respect to rhotor

        Input: -m or --miller
        Returns miller_r to rhotor relationship
        """)

    #Gradient Calculator with respect to rhotor#
    if (args.rho):

        print('Location   = {}'.format(rho_tr))

        print('omne_tr    = {}'.format(omne_tr[0]))
        print('omni_tr    = {}'.format(omni_tr[0]))
        print('omnb_tr    = {}'.format(omnb_tr[0]))
        print('omte_tr    = {}'.format(omte_tr[0]))
        print('omti_tr    = {}'.format(omti_tr[0]))
        print('omtb_tr    = {}'.format(omtb_tr[0]))
        print('omnc_tr    = {}'.format(omnc_tr[0]))

        print('nref_tr    = {}'.format(nref_tr))
        print('tref_tr    = {}'.format(tref_tr))
        print('nc_tr      = {}'.format(nc_tr))
        print('ne_tr      = {}'.format(ne_tr))
        print('ni_tr      = {}'.format(ni_tr))
        print('te_tr      = {}'.format(te_tr))
        print('ti_tr      = {}'.format(ti_tr))
        print('Tref       = {}'.format(tref_tr))

        print('me_tr      = {}'.format(me_tr))
        print('mD_tr      = {}'.format(mD_tr))
        print('rho_star   = {}'.format(rho_star))

        print('omeg_tr    = {}'.format(omeg_tr*1000))
        import extract_miller_from_eqdsk

    #Miller Geometery Relationship to rhotor#
    if (args.miller):
        mill = miller.findmiller(efitfpath)
        rad = mill[0]
        setparam = {'nrhomesh':'rhotor'}

        efitprof = cheasefiles.read_eqdsk(eqdskfpath=efitfpath,setParam=setparam)
        specprof = cheasefiles.read_profiles(profilesfpath=proffpath,setParam=setparam,eqdsk=efitfpath)
        rhopsi = np.array(specprof['rhopsi'])
        rhotor = np.array(specprof['rhotor'])

        plt.figure()
        plt.plot(rad,rhotor, linewidth=3.0)
        plt.ylabel('r/a', fontsize=16)
        plt.xlabel('$\\rho_{tor}$', fontsize=16)
        plt.grid()
        plt.show()

    #rho_star Calculation and Scaling of Profiles
    if (args.rho_star):
        iter_eqdsk = '/home/calderhaubrich/Workspace/GENESCRIPT/Cheasefiles/eqdsk.DT_15MA_Q10'
        iterpfile = '/home/calderhaubrich/Workspace/GENESCRIPT/Cheasefiles/ITER_DT_15MA_Q10_peaked_n.txt'

        def ReadPfile(filepath):
            f     = open(filepath,'r')
            pfile = f.readlines()
            del pfile[0:18]
            rawdata = []
            for i in range(len(pfile)):
                row = pfile[i].split()
                rawdata.append(row)
            profile_chunks = np.array(rawdata).T.tolist()
            for each in range(len(profile_chunks)):
                del profile_chunks[each][0]
                for num in range(len(profile_chunks[each])):
                    profile_chunks[each][num] = float(profile_chunks[each][num])
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
            plt.show()

        iter_case = millergeometryfunction.finder(iter_eqdsk, rho_tr)
        iter_Bref = iter_case[8]
        iter_Lref = iter_case[0]
        drho_tor_dr_iter = iter_case[2]
        radial = iter_case[7]

        #ITER Gradient Calculation
        iterprofiles = ReadPfile(iterpfile)
        K_boltz = 8.617 * 10 ** (-5) #Boltzmann Constant in eV/K

        T_i      = iterprofiles[10] #in kev
        n_W      = iterprofiles[15]
        rho_iter = iterprofiles[50]
        n_e      = iterprofiles[3]
        T_e      = iterprofiles[2]
        n_D      = iterprofiles[4]
        Vr       = iterprofiles[-1]

        T_i_spl = InterpolatedUnivariateSpline(rho_iter, T_i, k=interpol_order)
        n_e_spl = InterpolatedUnivariateSpline(rho_iter, n_e, k=interpol_order)
        T_e_spl = InterpolatedUnivariateSpline(rho_iter, T_e, k=interpol_order)
        n_D_spl = InterpolatedUnivariateSpline(rho_iter, n_D, k=interpol_order)
        Vrad_sec_spl = InterpolatedUnivariateSpline(rho_iter, Vr, k=interpol_order)

        dT_idrho_spl = T_i_spl.derivative()
        dn_edrho_spl = n_e_spl.derivative()
        dT_edrho_spl = T_e_spl.derivative()
        dn_Ddrho_spl = n_D_spl.derivative()

        t_ref = T_i_spl(rho_iter)
        Vr_loc = Vrad_sec_spl(rho_tr)
        nref_tr = n_e_spl(rho_tr)     # Electron density as refernce density at target location
        tref_tr = T_i_spl(rho_tr)
        ne_iter= n_e_spl(rho_tr)/nref_tr
        ni_iter= n_D_spl(rho_tr)/nref_tr
        #vtor_tr= vtor_spl(rho_tr)   # toroidal rotation velocity m/s at reference location
        #omeg_tr= omeg_spl(rho_tr)   # toroidal rotation frequency in rad/sec at reference location
        Vrad_sec = Vr_loc*1000/(iter_Lref+(radial*iter_Lref))
        nw_tr = (ne_iter - ni_iter)/37                                  # Tungsten Density with normalized mass

        T_e_tr= T_e_spl(rho_tr)/tref_tr
        T_i_tr= T_i_spl(rho_tr)/tref_tr

        omn_e_tr= dn_edrho_spl(rho_tr)*[-iter_Lref/n_e_spl(rho_tr)]*drho_tor_dr_iter #omne_tr = -(Lref/n)(dn/drho)
        omT_e_tr= dT_edrho_spl(rho_tr)*[-iter_Lref/T_e_spl(rho_tr)]*drho_tor_dr_iter
        omT_i_tr= dT_idrho_spl(rho_tr)*[-iter_Lref/T_i_spl(rho_tr)]*drho_tor_dr_iter
        omn_D_tr= dn_Ddrho_spl(rho_tr)*[-iter_Lref/n_D_spl(rho_tr)]*drho_tor_dr_iter
        omn_W_tr =(ne_iter*omn_e_tr -ni_iter*omn_D_tr)/(37*nw_tr)

        print('Location   = {}'.format(rho_tr))

        print('omne_tr    = {}'.format(omn_e_tr[0]))
        print('omni_tr    = {}'.format(omn_D_tr[0]))
        print('omnW_tr    = {}'.format(omn_W_tr[0]))
        print('omte_tr    = {}'.format(omT_e_tr[0]))
        print('omti_tr    = {}'.format(omT_i_tr[0]))
        print('nW_tr      = {}'.format(nw_tr))
        print('ni_tr      = {}'.format(ni_iter))
        print('ne_tr      = {}'.format(ne_iter))
        print('te_tr      = {}'.format(T_e_tr))
        print('ti_tr      = {}'.format(T_i_tr))
        print('Tref       = {}'.format(tref_tr))
        print('nref_tr    = {}'.format(nref_tr))
        print('omegatorref= {}'.format(Vrad_sec))

        #Gyroradius-to-machine-size ratio at reference location
        #rho_star = (m_ref*np.sqrt(tref_tr/m_ref))/(Bref*R0)
        #rho_star_tr = (m_ref*np.sqrt(t_ref/m_ref))/(Bref_tr*Lref_tr)

        #ITER Impurity Flux Calculation and Diffusion Coefficent#

        charge = []
        for T in range(len(T_i)):
            state = ionize.ionization(T_i[T])
            charge.append(state)
        #print(max(set(charge), key=charge.count))
        #print(sum(charge)/len(charge))
        #print(T_i)
        D_W = [] #Tungsten Impurity Bohm Diffusion Coefficent
        for tcharge in range(len(T_i)):
            T_i[tcharge] = T_i[tcharge] * 1000
            rati = (T_i[tcharge]/charge[tcharge])
            Diff = (1/16)*(K_boltz/abs(iter_Bref))*(rati)
            D_W.append(Diff)

        plt.figure()
        plt.plot(rho_iter,D_W, linewidth=3.0)
        plt.xlabel("R/$L_{n}$", fontsize=16)
        plt.ylabel('$D_{W}$(m^2/s)', fontsize=16)
        plt.tight_layout
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(rho_iter,T_i, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("T$_i$(eV)", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        fig,ax1 = plt.subplots()
        ax2=ax1.twinx()
        ax1.plot(rho_iter,T_i, 'g-')
        ax2.plot(rho_iter, charge, 'b-')
        ax1.set_xlabel("$\\rho_{tor}$")
        ax1.set_ylabel("T$_i$(eV)",color='g')
        ax2.set_ylabel("Z$_W$",color='b')
        plt.show()

        plt.figure()
        plt.plot(T_i,charge, linewidth=3.0)
        plt.xscale('log')
        plt.xlabel("T$_i$(eV)", fontsize=16)
        plt.ylabel("Z$_W$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        plt.figure()
        plt.plot(rho_iter, n_W, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("n$_W (m^{-3})$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

    #Species Plots with Respect to rhotor#
    if (args.all):
        # Ion temperature plot
        plt.figure()
        plt.plot(rhotor,Ti, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("T$_i$(eV)", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # Fast Ion temperature plot
        plt.figure()
        plt.plot(rhotor,Tb, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("T$_b$(eV)", fontsize=16)
        plt.tight_layout()
        plt.tight_layout()
        plt.grid()
        plt.show()


        # Electron temperature plot
        plt.figure()
        plt.plot(rhotor,Te, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("T$_e$(eV)", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()


        # Electron and ion temperature profiles plot
        plt.figure()
        plt.plot(rhotor,Te, linewidth=3.0)
        plt.plot(rhotor,Ti, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("T (eV)", fontsize=16)
        plt.legend(["e", "i"], prop={"size":16})
        plt.tight_layout()
        plt.grid()
        plt.show()


        # electron density plot
        plt.figure()
        plt.plot(rhotor,ne, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("n$_e (m^{-3})$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # Ion density plot
        plt.figure()
        plt.plot(rhotor,ni, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("n$_i (m^{-3})$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # NBI fast Ion density plot
        plt.figure()
        plt.plot(rhotor,nb, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$",fontsize=16 )
        plt.ylabel("n$_b (m^{-3})$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # NBI fast Ion density plot
        plt.figure()
        plt.plot(rhotor,Pb, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("P$_b (m^{-3})$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # Pressure profiles plot
        plt.figure()
        plt.plot(rhotor,Ptot, linewidth=3.0)
        plt.plot(rhotor,Pb, linewidth=3.0)
        plt.plot(rhotor,Ptot-Pb, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("P(Pa)", fontsize=16)
        plt.legend(["Total", "FI", "Th"], prop={"size":16})
        plt.tight_layout()
        plt.grid()
        plt.show()


        # Carbon density plot
        plt.figure()
        plt.plot(rhotor,nc, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("n$_C (m^{-3})$",  fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # Toroidal rotataion velocity
        plt.figure()
        plt.plot(rhotor,Vtor, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("V$_{tor}$ (m/s)",fontsize=16 )
        plt.tight_layout()
        plt.grid()
        plt.show()

        #Gyroradius-to-machine-size ratio plot
        plt.figure()
        plt.plot(rhotor,rho_star_tr, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("$\\rho_*$", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # Safety factor plot
        plt.figure()
        plt.plot(rhotor,q, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("Safety factor (q)", fontsize=16)
        plt.tight_layout()
        plt.grid()
        plt.show()

        # electron and ion density plot
        plt.figure()
        plt.plot(rhotor,ne, linewidth=3.0)
        plt.plot(rhotor,ni, linewidth=3.0)
        plt.xlabel("$\\rho_{tor}$", fontsize=16)
        plt.ylabel("n $(m^{-3})$", fontsize=16)
        plt.legend(["e", "i"], prop={"size":16})
        plt.tight_layout()
        plt.grid()
        plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])
