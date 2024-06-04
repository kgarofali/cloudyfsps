#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from past.utils import old_div
from builtins import range
from __future__ import (division, print_function, absolute_import, unicode_literals)

import sys
import numpy as np
import fsps
from cloudyfsps.ASCIItools import (writeASCII, compileASCII, checkCompiled, compiledExists)

# this code snippet goes through every step needed
# to integrate FSPS into Cloudy.
# This example uses stellar pops with a constant SFH
# as the input ionizing source.
# 1. Write an ascii file in Cloudy format with grid
#    of FSPS spectra in all available ages and
#    metallicities
# 2. Compile asii file into binary format required
#    for Cloudy use. Assumes $CLOUDY_EXE is set to
#    your /path/to/cloudy.exe
# 3. Writes Cloudy input files for a subset of grid
#    parameters.
# 4. Runs Cloudy on the *.in files
# 5. Formats the various output files

zsun = 0.02 # this is solar metallicity for the MIST isochrones
# in BPASS models zsun = 0.02
exec_write_ascii = True

# Function to write the ascii file.
# This is where you set the properties of the
# ionizing spectrum (SSP/CSFH, IMF, FBHB, etc)

def mist_ascii(fileout, **kwargs):
    # change these parameters to modify the ionizing source grid
    # default mode is to produce an ascii grid in age and Z,
    # though different variables and more dimensions are possible.

    # or sfh = 0 for SSP (see below, i guess this would be simple burst?)
    # can set dust_type=2 for calzetti or 5 for SMC
    # imf_type=0 for salpeter

    sp_dict = dict(zcontinuous=1,
                   imf_type=2,
                   sfh=0,
                   const=0.0,
                   sf_start=0.0)
    sp = fsps.StellarPopulation(**sp_dict)
    # all ages and Zs
    all_ages = 10.**sp.log_age
    logZs = np.log10(old_div(sp.zlegend,zsun))
    lam = sp.wavelengths
    # taking these as rounded to nearest decimal place closest values in SSP full grid for ages of 1, 3, 5, 7, 10, 15, and 20 Myr (but is it reasonable to assume we still have a cloud at 20 Myr?)
    # where XRB will be necessarily zero at 1 & 3 Myr, same value from earliest times in Fragos grid (~ 7 Myr) at 5 and 7 Myr, then interpolated values at 10, 15, and 20 Myr
    ages_grid = np.array([1.0,3.2,5,7.9,10,15.8,20])
    ages_idx = []
    for age in ages_grid:
        rounded_all_ages = np.array([round(step, 1) for step in (all_ages/1e+6)])
        match_idx, = np.where(age == rounded_all_ages)
        ages_idx.append(int(match_idx))
    # just run this at 0.1 Zsun and solar for consistency with values from XRB scalings.
    logZs_idx = []
    all_fluxs = []
    for logZ in logZs:
        if (logZ == -1) | (logZ == 0):
            match_idx, = np.where(logZ == logZs)
            logZs_idx.append(int(match_idx))
            sp.params['logzsol'] = logZ
            all_fluxs.append(sp.get_spectrum()[1]) #lsun per hz
    modpars = [(age, logZ) for age in all_ages[ages_idx] for logZ in logZs[logZs_idx]]
    nmod = len(modpars)
    # flatten flux for writing
    flat_flux = np.array([all_fluxs[j][i] for i in ages_idx for j in range(len(logZs[logZs_idx]))])

    # figure out stellar mass scaling at each grid point
    # arr in Qh_hat (as defined in Byler+2017 S2.1.4) will be 0.1 Zsun 1,3,5,7,10,15,20 Myr, followed by Zsun 1,3,5,7,10,15,20 Myr.
    Qh_hat = np.array([calcQ(lam,all_fluxs[j][i]*3.839e+33,f_nu=True) for i in ages_idx for j in range(len(logZs[logZs_idx]))])
    # for our arr of logU from -4 --> -1 in bins of 0.5, for a fixed radius of 10**19 cm and a fixed density of 100 cm^-3, using eq. 2 in Byler+2017
    Uo_cloudy_input = np.logspace(-4,-1,7)
    Qh_cloudy_input = Uo_cloudy_input*4*np.pi*(10**19)**2*100*2.99e+10
    # Then we want to evaluate Qh_cloudy_input at each Qh_hat, so that is 7 instances of Qh_cloudy_input to evaluate for our age/metallicity grid of 14 points, or 98 independent grid points
    Qh_grid = []
    for i in range(len(Qh_cloudy_input)):
        arr_to_append = np.zeros(len(Qh_hat)+2)
        arr_to_append[0] = Uo_cloudy_input[i]
        arr_to_append[1] = Qh_cloudy_input[i]
        arr_to_append[2:] = Qh_hat
        Qh_grid.append(arr_to_append)
    np.savetxt('qgrid.txt', dat, delimiter=' ', header='logU,Qh,Qhat(0.1,1),Qhat( 0.1,3),Qhat(0.1,5),Qhat(0.1,7),Qhat(0.1,10),Qhat(0.1,15),Qhat(0.1,20),Qhat(1,1),Qhat(1,3),Qhat(1,5),Qhat(1,7),Qhat(1,10),Qhat(1,15),Qhat(1,20)')

    # do this part for proper XRB scaling (not astropy not loaded in the pyfsps enviro)
    # these file paths are prob not right, doing this now for illustrative purposes
    qtab = ascii.read('qgrid.txt', names=['logU','Qh','Qhat(0.1,1)','Qhat(0.1,3)','Qhat(0.1,5)','Qhat(0.1,7)','Qhat(0.1,10)','Qhat(0.1,15)','Qhat(0.1,20)','Qhat(1,1)','Qhat(1,3)','Qhat(1,5)','Qhat(1,7)','Qhat(1,10)','Qhat(1,15)','Qhat(1,20)'])
    tage, lxt_mstar_sol = np.loadtxt('/Users/kgarofal/research/heii/xrb_cloudy_grids/code/fragos_tot_zsol.txt', unpack=True)
    tage, lxt_mstar_sub = np.loadtxt('/Users/kgarofal/research/heii/xrb_cloudy_grids/code/fragos_tot_01zsol.txt', unpack=True)
    tage, lxt_mstar_sol_hmxb = np.loadtxt('/Users/kgarofal/research/heii/xrb_cloudy_grids/code/fragos_hmxb_zsol.txt', unpack=True)
    lxt_mstar_sub_hmxb = (lxt_mstar_sol*lxt_mstar_sub)*(lxt_mstar_sol_hmxb/lxt_mstar_sol)
    # from scipy import interpolate
    f_interp_sub_hmxb = interpolate.interp1d(np.log10(tage*1e+9),np.log10(lxt_mstar_sub_hmxb))
    f_interp_sol_hmxb = interpolate.interp1d(np.log10(tage*1e+9),np.log10(lxt_mstar_sol_hmxb))
    # arr in lx_mstar_hmxb will be 0.1 Zsun 1,3,5,7,10,15,20 Myr, followed by Zsun 1,3,5,7,10,15,20 Myr as for Qh_hat
    lxt_mstar_hmxb_grid = np.zeros(len(Qh_hat))
    # first fill out subsolar and solar values (indexed by 7)
    for i in range(len(ages_grid)):
        if ages_grid[i] < 5.:
            lxt_mstar_hmxb_grid[i] = 0.
            lxt_mstar_hmxb_grid[i+7] = 0.
        elif (ages_grid[i] >= 5.) & (ages_grid[i] < 7):
            lxt_mstar_hmxb_grid[i] = 10**f_interp_sub_hmxb(np.log10(ages_grid[i+1]*1e+6))
            lxt_mstar_hmxb_grid[i+7] = 10**f_interp_sol_hmxb(np.log10(ages_grid[i+1]*1e+6))
        else:
            lxt_mstar_hmxb_grid[i] = 10**f_interp_sub_hmxb(np.log10(ages_grid[i]*1e+6))
            lxt_mstar_hmxb_grid[i+7] = 10**f_interp_sol_hmxb(np.log10(ages_grid[i]*1e+6))

    # the answers i get from this seem rather low in luminosity, because we aren't forming much stellar mass for a given Uo/Qh. maybe that's accurate, but if we choose this route by definiton we won't get much contribution i would imagine...
    lxbol_tscal_hmxb = []
    lxband_empscal_hmxb = []
    sfrs = []
    # average lx/sfr for 12+logOH = 8.2 and 8.6 from Lehmer+2021
    lx_sfr_avg = np.zeros(len(lxt_mstar_hmxb_grid))
    lx_sfr_avg[:7] = 10**39.94
    lx_sfr_avg[7:] = 10**39.64
    for i in range(len(lxt_mstar_hmxb_grid)):
        # skip first two colums indexing table because these just record Uo and Qh. first term in this the total stellar mass for a given Uo/Qh, and second term is the lx per unit mass at that age and metallicity (hopefully properly indexed)
        lxbol_scal_per_Uo = ((qtab['Qh']/qtab.columns[i+2])/1e+10)*lxt_mstar_hmxb_grid[i]
        lxbol_tscal_hmxb.append(lxbol_scal_per_Uo)
        sfrs.append((qtab['Qh']/qtab.columns[i+2])/(concatenate((ages_grid, ages_grid))[i]*1e+6))
        # mult by average lx/sfr value
        lxband_empscal_hmxb.append(((qtab['Qh']/qtab.columns[i+2])/(concatenate((ages_grid, ages_grid))[i]*1e+6))*(lx_sfr_avg[i]))

    ## hanging on to some code for plotting functionality
    ## see here: https://stackoverflow.com/questions/26545897/drawing-a-colorbar-aside-a-line-plot-using-matplotlib
    fig, ax = plt.subplots(figsize=(8,8))
    cmap_params = np.log10(Uo_cloudy_input)
    norm = matplotlib.colors.Normalize(vmin=np.min(cmap_params),vmax=np.max(cmap_params))
    c_m = matplotlib.cm.plasma
    s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
    s_m.set_array([])
    for i in range(len(qtab)):
        # with theory based scaling based on stellar mass formed
        ax.plot(np.log10(ages_grid*1e+6),np.log10(((qtab[i]['Qh']/np.array(list(qtab[i]['Qhat(1,1)','Qhat(1,3)','Qhat(1,5)','Qhat(1,7)','Qhat(1,10)','Qhat(1,15)','Qhat(1,20)'])))/1e+10)*lxt_mstar_hmxb_grid[7:]),'--', lw=2, color=s_m.to_rgba(cmap_params[i]))
        #with empirical scaling based on avg sfr
        ax.plot(np.log10(ages_grid*1e+6),np.log10(((qtab[i]['Qh']/np.array(list(qtab[i]['Qhat(1,1)','Qhat(1,3)','Qhat(1,5)','Qhat(1,7)','Qhat(1,10)','Qhat(1,15)','Qhat(1,20)'])))/(ages_grid*1e+6))*lx_sfr_avg[7:]),'-', lw=2, color=s_m.to_rgba(cmap_params[i]))
    cbar = plt.colorbar(s_m)
    cbar.set_label('logU')
    ax.set_xlabel('log(Age)')
    ax.set_ylabel('log(Lx)')
    plt.show()

    # this function is flexible, ndim can be 3/4/n.
    # in this example, however, ndim is 2 (age, logz).
    # note for me this writes to ~/software/cloudy/c17.02/data
    writeASCII(fileout, lam, flat_flux, modpars,
               nx=len(lam), ndim=2, npar=2, nmod=nmod)
    return
#---------------------------------------------------------------------
# ASCII FILE: WRITE AND COMPILE
#---------------------------------------------------------------------
# assumes you have $CLOUDY_EXE and $CLOUDY_DATA_PATH set as sys vars.

# name of ascii file
ascii_file = 'FSPS_MIST_SSP.ascii'

# the ascii file takes a while to generate, so if an already-compiled
# version exists, the code will not overwrite it.

compiled_ascii = '{}.mod'.format(ascii_file.split('.')[0])
if exec_write_ascii:
    print("Executing write ascii sequence...")
    if not compiledExists(ascii_file):
        print("No compiled model exists...Writing.")
        mist_ascii(ascii_file)
        print("Compiling {} with Cloudy".format(ascii_file))
        compileASCII(ascii_file)
        print("Checking to see if compilation was successful...")
        if checkCompiled(ascii_file):
            print("Your model {} is ready to run.".format(compiled_ascii))
        else:
            sys.exit()
    else:
        print("{} already exists.".format(compiled_ascii))
