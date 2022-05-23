# -*- coding: utf-8 -*-

import cheasefiles
import millergeometryfunction
import rho_tor_to_miller_r as miller
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline



plt.rcParams.update({'font.size':14})

proffpath = '/home/calderhaubrich/GENE/DIIID_163241_profiles/p163241.03500_e5059'

def main(args):
    parser = argparse.ArgumentParser(description='Calculates gradients at given rhotor position.')
    parser.add_argument('file', help='input some gfile')
    parser.add_argument('location', help='rhotor location')
    parser.add_argument('-r','--rho', default=False, action='store_true', help='Finds rho gradients')
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
        geomdata = millergeometryfunction.finder(efitfpath,rho_tr)
        R0 = geomdata[0]
        a  = geomdata[1]
        drho_tor_dr = geomdata[2]

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

        dnedrho = dnedrho_spl(rhotor)
        dnidrho = dnidrho_spl(rhotor)
        dnbdrho = dnbdrho_spl(rhotor)
        dtedrho = dtedrho_spl(rhotor)
        dtidrho = dtidrho_spl(rhotor)
        dtbdrho = dtbdrho_spl(rhotor)
        dpbdrho = dpbdrho_spl(rhotor)
        dvtordrho = dvtordrho_spl(rhotor)
        #domegarho = domegdrho_spl(rhotor)

        #dvpolrho = dvpoldrho_spl(rhotor)


        # Normalized paramters at target toroidal radial location
        Z0= 0.0031727/R0
        mD= 3.343583719e-27         # Mass of deuterium in kg
        m_ref= mD                   # Deuteriam mass as reference mass
        mD_tr= mD/m_ref            # Normalized deuterium mass
        me= 9.10938356e-31          # Mass of electron in SI Units (kg)
        me_tr= me/m_ref             # Normalized electron mass
        #n_ref= ne_spl(rhotor)      # reference temparture as function of rhotor
        #t_ref= ti_spl(rhotor) 
        nref_tr= ne_spl(rho_tr)     # Electron density as refernce density at target location
        tref_tr= ti_spl(rho_tr)     # main ion tempeartrure as reference temperature at target location

        vtor_tr= vtor_spl(rho_tr)     # toroidal rotation velocity m/s at reference location
        omeg_tr= omeg_spl(rho_tr)    # toroidal rotation frequency in rad/sec at reference location

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

        #Calculate Gradients with respect to minor radius a
        omne = dnedrho_spl(rho_tr)*[-a/ne_spl(rho_tr)]*drho_tor_dr
        omni = dnidrho_spl(rho_tr)*[-a/ni_spl(rho_tr)]*drho_tor_dr
        omnb = dnbdrho_spl(rho_tr)*[-a/nb_spl(rho_tr)]*drho_tor_dr
        omte = dtedrho_spl(rho_tr)*[-a/te_spl(rho_tr)]*drho_tor_dr
        omti = dtidrho_spl(rho_tr)*[-a/ti_spl(rho_tr)]*drho_tor_dr
        omtb = dtbdrho_spl(rho_tr)*[-a/tb_spl(rho_tr)]*drho_tor_dr

        print('Location   = {}'.format(rho_tr))

        print('omne_tr    = {}'.format(omne_tr))
        print('omni_tr    = {}'.format(omni_tr))
        print('omnb_tr    = {}'.format(omnb_tr))
        print('omte_tr    = {}'.format(omte_tr))
        print('omti_tr    = {}'.format(omti_tr))
        print('omtb_tr    = {}'.format(omtb_tr))
        print('omnc_tr    = {}'.format(omnc_tr))

        print('nref_tr    = {}'.format(nref_tr))
        print('tref_tr    = {}'.format(tref_tr))
        print('nc_tr      = {}'.format(nc_tr))
        print('ne_tr      = {}'.format(ne_tr))
        print('ni_tr      = {}'.format(ni_tr))
        print('te_tr      = {}'.format(te_tr))
        print('ti_tr      = {}'.format(ti_tr))

        print('me_tr      = {}'.format(me_tr))
        print('mD_tr      = {}'.format(mD_tr))

        print('omeg_tr    = {}'.format(omeg_tr*1000))

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
