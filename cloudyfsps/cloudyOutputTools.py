#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from builtins import str as newstr
from builtins import range
#__all__ = ["format_output"]

import numpy as np
import subprocess
import pkg_resources
from .generalTools import air_to_vac
from scipy.interpolate import interp1d
import pandas as pd
###
# ***.lin: [cloudy_ID, flux]
# ***.lineflux: [sorted_vac_wl, flux]
# ***.out_lines: [sorted_vac_wl, flux]
###
# ***.outwcont: [wl, attenuated_incident, diffuse_continuum]
# ***.inicont: [wl, incident_flux]
# ***.contflux: [wl, incid_out, atten_out, diffuse_out]
# ***.out_cont: [ang, diffuse_out]
###
def formatCloudyOutput(dir_, model_prefix, modnum, modpars, use_extended_lines=True, write_line_lum=True, **kwargs):
    '''
    for formatting the output of a single cloudy job
    '''
    # model information
    logZ, age, logU, logR, logQ, nH = modpars[0:6]
    if logZ > 0.2:
        print("WARNING WARNING WARNING")

    dist_fact = 4.0*np.pi*(10.0**logR)**2.0 # cm**2
    lsun = 3.839e33 # erg/s
    c = 2.9979e18 #ang/s

    oldfile = "{}{}{}.lin".format(dir_, model_prefix, str(modnum))
    newfile = "{}{}{}.lineflux".format(dir_, model_prefix, str(modnum))
    print_file = "{}{}{}.out_lines".format(dir_, model_prefix, str(modnum))
    # read cloudy output
    dat = np.genfromtxt(oldfile, skip_header=2, delimiter="\t",
                        dtype="S20,f8")
    #line_names = [d[0] for d in dat]
    # read in fluxes and line namesfrom cloudy .lin file
    # strip the line names to be just the wavelength, and convert all to angstroms
    datflu = np.array([d[1] for d in dat])
    dat_line_names = np.array([d[0] for d in dat])
    datwl = []
    for line_name in dat_line_names:
        ln_str = str(line_name).split()[-2]
        if ln_str[-1] == 'm':
            ln_wl = float(ln_str[:-1])*1e+4
        elif ln_str[-1] == 'A':
            ln_wl = float(ln_str[:-1])
        datwl.append(ln_wl)
    datwl = np.array(datwl)
    # non-ordered wavelengths
    if use_extended_lines:
        wavfile = pkg_resources.resource_filename(__name__,
                                                  "data/refLinesEXT.dat")
    else:
        wavfile = pkg_resources.resource_filename(__name__,
                                                  "data/refLines.dat")
    wdat = np.genfromtxt(wavfile, delimiter=',', dtype=None)
    wl = np.array([dat[0] for dat in wdat])
    # sort them by wavelength
    sinds = np.argsort(wl)

    output = np.column_stack((wl[sinds], datflu[sinds]))
    np.savetxt(newfile, output, fmt=str("%4.6e"))
    # print lines to ***.out_lines
    # lines are originall saved from cloudy in units of erg/s, so "conv" option gives them in units of solar lums per Q, otherwise, prints them as absolute luminosities
    line_wav = wl[sinds]
    if write_line_lum:
        conv = 1.0
    else:
        conv = 1./lsun/(10.**logQ)
    line_flu = datflu[sinds]*conv
    print_output = np.column_stack((line_wav, line_flu))
    np.savetxt(print_file, print_output, fmt=(str("%.6e"),str("%.6e")))
    # print to file
    print("Lines were printed to file {}".format(print_file))
    ########
    ### continuum
    ########
    outcontfl = "{}{}{}.outwcont".format(dir_, model_prefix, modnum)
    incontfl = "{}{}{}.inicont".format(dir_, model_prefix, modnum)
    print_file2 = "{}{}{}.contflux".format(dir_, model_prefix, modnum)
    print_file = "{}{}{}.out_cont".format(dir_, model_prefix, modnum)
    # lam, atten_inc, diff_cont, diff_line, sum
    cont_df = pd.read_csv(outcontfl, skiprows=1, delimiter='\t')
    cont_data = cont_df.values
    # BIG difference from Nell's code, cont (and other emission properties from e.g., save continuum) are in nuLnu, that is Hz * (erg/s/Hz), or erg/s
    # Thus, we will write things out in terms of absolute luminosities, so if we want to get things in terms of a flux later need to divide by dist_factor, per unit wavelength divide by angstrom, per solar lum divide by lsun, etc.
    atten_0, diffuse_0, cloud_tot_0, lines_0  = cont_data[:,2], cont_data[:,3], cont_data[:,4], cont_data[:,8]
    cont_0 = cont_data[:,3] - cont_data[:,8]
    ang_0 = cont_data[:,0]
    # reverse arrays
    atten_in, diffuse_in, cloud_tot_in, lines_in, cont_in = atten_0[::-1], diffuse_0[::-1], cloud_tot_0[::-1], lines_0[::-1], cont_0[::-1]
    ang = ang_0[::-1]
    ang_v = air_to_vac(ang)
    # interpolate
    lamfile = pkg_resources.resource_filename(__name__, "data/FSPSlam.dat")
    fsps_lam = np.genfromtxt(lamfile)
    nu = c/fsps_lam
    atten_y = interp1d(ang_v, atten_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    diffuse_y = interp1d(ang_v, diffuse_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    cloud_tot_y = interp1d(ang_v, cloud_tot_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    cont_y = interp1d(ang_v, cont_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    ##
    # diffuse continuum in Lsun/Qh/Hz for use with fsps and prospector modelling
    diffuse_out = (diffuse_y) / (nu*(10.**logQ)*lsun)
    cont_out = (cont_y) / (nu*(10.**logQ)*lsun)
    ##
    inidata = np.genfromtxt(incontfl, skip_header=1)
    angini_0 = inidata[:,0]
    angini = angini_0[::-1]
    angini_v = air_to_vac(angini)
    incid_0 = inidata[:,1]
    incid_in = incid_0[::-1]
    incid_y = interp1d(angini_v, incid_in, fill_value=0.0, bounds_error=False)(fsps_lam)
    # F_nu / (nu=c/lambda) per solar lum
    f = open(print_file2, "w")
    f.write("# lam (ang) incid (erg/s) attenuated_incid (erg/s) neb_tot neb_cont (erg/s)\n")
    for i in range(len(fsps_lam)):
        printstring = "{0:.6e} {1:.6e} {2:.6e} {3:.6e} {4:.6e}\n".format(fsps_lam[i], incid_y[i], atten_y[i], diffuse_y[i], cont_y[i])
        f.write(printstring)
    f.close()
    print("The full continuum was printed to file {}".format(print_file2))
    #####
    f = open(print_file, "w")
    f.write("# lam (ang) diffuse_cont (lsun/hz/Q)\n")
    for i in range(len(fsps_lam)):
        printstring = "{0:.6e} {1:.6e} {2:.6e}\n".format(fsps_lam[i], diffuse_out[i], cont_out[i])
        f.write(printstring)
    f.close()
    print("The diffuse continuum was printed to file {}".format(print_file))
    return

def formatAllOutput(dir_, mod_prefix, use_extended_lines=True, write_line_lum=True):
    '''
    for formatting output after running a batch of cloudy jobs
    '''
    data = np.genfromtxt(dir_+mod_prefix+".pars")
    def get_pars(modnum):
        return data[np.int(modnum)-1, 1:]
    for modnum in data[:,0]:
        mnum = np.int(modnum)
        formatCloudyOutput(dir_, mod_prefix, mnum, get_pars(mnum), use_extended_lines=use_extended_lines, write_line_lum=write_line_lum)
    return
