#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cheasefiles
import millergeometryfunction
import numpy as np
import sys
import argparse
from scipy.interpolate import InterpolatedUnivariateSpline

def main(args):
    parser = argparse.ArgumentParser(description='Creates Parameters file for GENE')
    parser.add_argument('gfile', help='input some gfile')
    parser.add_argument('pfile', help='input some pfile')
    parser.add_argument('location', help='rhotor location')
    parser.add_argument('-i','--info', default=False, action='store_true', help='Provides info on code')
    args = parser.parse_args(args)

    efitfpath = args.gfile
    proffpath = args.pfile
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
    R0     = geomdata[0]
    a      = geomdata[1]
    drho_tor_dr = geomdata[2]
    q0     = geomdata[5]
    shat   = geomdata[6]
    trpeps = geomdata[7]
    Bref   = geomdata[8]
    amhd   = geomdata[9]
    drR    = geomdata[10]
    drZ    = geomdata[11]
    kappa  = geomdata[12]
    s_kappa= geomdata[13]
    delta  = geomdata[14]
    s_delta= geomdata[15]
    zeta   = geomdata[16]
    s_zeta = geomdata[17]
    minor_r= geomdata[18]

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

    #User inputs for parameters#
    #parallelization#
    n_procs = input("Input values for n_procs(s v w x y z sim): ")
    if n_procs == 'exit':
        quit()
    n_procs_list = list(map(int, n_procs.split()))
    while len(n_procs_list) != 7:
        print("Error: n_procs requires 7 entries.")
        n_procs = input("Input values for n_procs as list: ")
        n_procs_list = list(map(int, n_procs.split()))
        if n_procs == 'exit':
            quit()
    n_procs_s   = int(n_procs_list[0])
    n_procs_v   = int(n_procs_list[1])
    n_procs_w   = int(n_procs_list[2])
    n_procs_x   = int(n_procs_list[3])
    n_procs_y   = int(n_procs_list[4])
    n_procs_z   = int(n_procs_list[5])
    n_procs_sim = int(n_procs_list[6])

    #box#
    box1 = input("Input values for n_spec nx0 nky0 nz0 nv0 nw0: ")
    if box1 == 'exit':
        quit()
    box1_list = list(map(int, box1.split()))
    while len(box1_list) != 6:
        print("Error: requires 6 entries.")
        box1 = input("Input values for n_spec nx0 nky0 nz0 nv0 nw0: ")
        box1_list = list(map(int, box1.split()))
        if box1 == 'exit':
            quit()
    n_spec = int(box1_list[0])
    nx0    = int(box1_list[1])
    nky0   = int(box1_list[2])
    nz0    = int(box1_list[3])
    nv0    = int(box1_list[4])
    nw0    = int(box1_list[5])

    box2 = input("Input values for kymin lv lw lx: ")
    if box2 == 'exit':
        quit()
    box2_list = list(map(float, box2.split()))
    while len(box2_list) != 4:
        print("Error: requires 4 entries.")
        box2 = input("Input values for kymin lv lw lx: ")
        box2_list = list(map(float, box2.split()))
        if box2 == 'exit':
            quit()
    kymin = box2_list[0]
    lv    = box2_list[1]
    lw    = box2_list[2]
    lx    = box2_list[3]

    #istep#
    istep = input("Input values for istep(field mom nrg omega vsp schpt energy): ")
    if istep == 'exit':
        quit()
    istep_list = list(map(int, istep.split()))
    while len(istep_list) != 7:
        print("Error: requires 7 entries.")
        istep = input("Input values for istep(field mom nrg omega vsp schpt energy): ")
        istep_list = list(map(int, istep.split()))
        if istep == 'exit':
            quit()
    istep_field  = int(istep_list[0])
    istep_mom    = int(istep_list[1])
    istep_nrg    = int(istep_list[2])
    istep_omega  = int(istep_list[3])
    istep_vsp    = int(istep_list[4])
    istep_schpt  = int(istep_list[5])
    istep_energy = int(istep_list[6])

    #perf_vec#
    perf_vec = input("Input perf_vec: ")
    if perf_vec == 'exit':
        quit()
    
    #nblocks#
    nblocks = input("Input nblocks: ")
    if nblocks == 'exit':
        quit()
    nblocks = int(nblocks)

    #time#
    timelim = input("Input timelim: ")
    if timelim == 'exit':
        quit()
    timelim = int(timelim)
    ntimesteps = input("Input ntimesteps: ")
    if ntimesteps == 'exit':
        quit()
    ntimesteps = int(ntimesteps)

    #parameters file#
    le = "\n"
    print('     ')
    print('&parallelization')
    print('n_procs_s   =   {}'.format(n_procs_s))
    print('n_procs_v   =   {}'.format(n_procs_v))
    print('n_procs_w   =   {}'.format(n_procs_w))
    print('n_procs_x   =   {}'.format(n_procs_x))
    print('n_procs_y   =   {}'.format(n_procs_y))
    print('n_procs_z   =   {}'.format(n_procs_z))
    print('n_procs_sim =   {}'.format(n_procs_sim))
    print('/'+le)

    print('&box')
    print('n_spec =   {}'.format(n_spec))
    print('nx0    =   {}'.format(nx0))
    print('nky0   =   {}'.format(nky0))
    print('nz0    =   {}'.format(nz0))
    print('nv0    =   {}'.format(nv0))
    print('nw0    =   {}'.format(nw0)+le)

    print('x0     =   {}'.format(rho_tr))
    print('kymin  =   {}'.format(kymin))
    print('lv     =   {}'.format(lv))
    print('lw     =   {}'.format(lw))
    print('lx     =   {}'.format(lx)+le)

    print('adapt_lx = T')
    print('adapt_ly = T')
    print('n0_global = -1111')
    print('/'+le)

    print('&in_out')
    print("diagdir = './'"+le)

    print('read_checkpoint  = F')
    print('write_checkpoint = T'+le)

    print('istep_field  =   {}'.format(istep_field))
    print('istep_mom    =   {}'.format(istep_mom))
    print('istep_nrg    =   {}'.format(istep_nrg))
    print('istep_omega  =   {}'.format(istep_omega))
    print('istep_vsp    =   {}'.format(istep_vsp))
    print('istep_schpt  =   {}'.format(istep_schpt))
    print('istep_energy =   {}'.format(istep_energy))
    print('/'+le)

    print('&general')
    print('nonlinear =   F')
    print("comp_type = 'IV'")
    print('perf_vec  =   {}'.format(perf_vec))
    print('nblocks   =   {}'.format(nblocks))
    print('x_local      = T')
    print('calc_dt      = .t.'+le)

    print('dt_max     =    0.000')
    print('dt_vlasov  =   0.7340E-03')
    print('ev_coll    =    2.3299')
    print('courant    =     1.25'+le)

    print('timelim    =   {}'.format(timelim))
    print('ntimesteps =   {}'.format(ntimesteps))
    print("omega_prec    = 1e-3")
    print("timescheme = 'RK4'")
    print('simtimelim =   0.1000E+05'+le)

    print("include_f0_contr = F"+le)

    print("init_cond = 'alm'"+le)

    print('beta       =   -1')
    print("debye2     =   -1")
    print("collision_op = 'landau'")
    print('coll       =   -1')
    print("coll_cons_model  = 'self_adj'"+le)

    print('hyp_z_with_dz_prefactor = F')
    print('hyp_z =    -1')
    print('hyp_x =      0.0')
    print('hyp_y =      0.0')
    print('hyp_v_with_dv_prefactor = F')
    print('hyp_v =   0.2000'+le)

    print('perf_tsteps =  -1'+le)

    print('/'+le)

    print('&external_contr')
    print('ExBrate    = 0'+le)

    print('pfsrate    = 0')
    print('Omega0_tor = -1111.0'+le)

    print('with_coriolis       = T')
    print('with_centrifugal    = T')
    print('with_comoving_other = T'+le)

    print('&geometry')
    print("magn_geometry = 'miller'"+le)

    print('rhostar = -1')

    print('q0       =   {}'.format(f"{q0:.5}"))
    print('shat     =   {}'.format(f"{shat:.5}"))
    print('trpeps   =   {}'.format(f"{trpeps:.5}"))
    print('amhd     =   {}'.format(f"{amhd:.5}"))
    print('drR      =   {}'.format(f"{drR:.5}"))
    print('drZ      =   {}'.format(f"{drZ:.5}"))
    print('kappa    =   {}'.format(f"{kappa:.5}"))
    print('s_kappa  =   {}'.format(f"{s_kappa:.5}"))
    print('delta    =   {}'.format(f"{delta:.5}"))
    print('s_delta  =   {}'.format(f"{s_delta:.5}"))
    print('zeta     =   {}'.format(f"{zeta:.5}"))
    print('s_zeta   =   {}'.format(f"{s_zeta:.5}"))
    print('minor_r  =   {}'.format(f"{minor_r:.5}"))
    print('major_R  =   1'+le)

    print('sign_Ip_CW =       1')
    print('sign_Bt_CW =       1')
    print('/'+le)

    print('&species')
    print("name   = 'ions'")
    print('omn    =  {}'.format(f"{omni_tr[0]:.9}"))
    print('omt    =  {}'.format(f"{omti_tr[0]:.9}")+le)

    print('mass   =  1.000')
    print('temp   =  1.000')
    print('dens   =  {}'.format(ni_tr))
    print('charge =  1')
    print('/'+le)

    print('&species')
    print("name   = 'electrons'")
    print('omn    =  {}'.format(f"{omne_tr[0]:.9}"))
    print('omt    =  {}'.format(f"{omte_tr[0]:.9}")+le)

    print('mass   =  {}'.format(me_tr))
    print('temp   =  {}'.format(te_tr))
    print('dens   =  1')
    print('charge = -1')
    print('/'+le)

    print('&species')
    print("name   = 'carbon'")
    print('omn    =  {}'.format(f"{omnc_tr[0]:.9}"))
    print('omt    =  {}'.format(f"{omti_tr[0]:.9}")+le)

    print('mass   =  6.000')
    print('temp   =  1.000')
    print('dens   =  {}'.format(nc_tr))
    print('charge =  6')
    print('/'+le)

    print('&units')
    print('Bref  =  {}'.format(f"{abs(Bref):.5}"))
    print('Tref  =  {}'.format(f"{tref_tr/1000:.5}"))
    print('nref  =  {}'.format(f"{nref_tr/1e+19:.5}"))
    print('Lref  =  {}'.format(f"{R0:.5}"))
    print('mref  =  1.99')
    print('omegatorref = {}'.format(omeg_tr*1000))
    print('/')

    if (args.info):
        print("""
Input gfile, pfile, and rho_tor location. This code creates and exports parameters file for GENE input.
        """)

if __name__ == "__main__":
    main(sys.argv[1:])
