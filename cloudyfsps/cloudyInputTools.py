#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import range
import os
import sys
import cloudyfsps
import numpy as np
import subprocess
import pkg_resources
from .generalTools import calcForLogQ
from .nebAbundTools import getNebAbunds

def cloudyInput(dir_, model_name, **kwargs):
    '''
    cloudyInput('./test/', 'ZAU115', logZ=-1.5, age=5.0e6, logU=-4.0)
    writes standard cloudy input to ./test/ZAU115.in
    defaults: 1Myr, logZ=-0.5, logU=-2.0, nH=100, r_inner=3 pc
    '''
    pars = {"age":1.0e6, #age in years
            "logZ": -0.5, #logZ/Zsol (-2.0 to 0.2)
            "gas_logZ":None,
            "logQ":47.0,
            "logU":-2.0, #log of ionization parameter
            "dens":100.0, # number density of hydrogen
            "r_inner":1.0, #inner radius of cloud
            "r_in_pc":False,
            "use_Q":True,
            "set_name":"dopita",
            "dust":True,
            "re_z":False,
            "cloudy_mod":"FSPS_SPS.mod",
            "efrac":-1.0,
            "extras":"",
            "extra_output":True,
            "to_file":True,
            "verbose":False,
            "par1":"age",
            "par1val":5.0e6,
            "par2":"logz",
            "par2val":0.0,
            "maxStellar":None,
            "use_extended_lines":False,
            "geometry":"sphere",
            "table_sed":False,
            "table_sed_frac": 0.0
            }
    for key, value in list(kwargs.items()):
        pars[key] = value
    # -----
    if pars["to_file"]:
        file_name = dir_+model_name+".in"
        f = open(file_name, "w")
    def this_print(s, eol=True):
        if s is None:
            print('"None" parameter not printed')
        else:
            to_print = s.strip()
            if pars["verbose"]:
                print(to_print)
            if pars["to_file"]:
                if eol: to_print += "\n"
                f.write(to_print)
    #-----------------------------------------------------
    if pars["gas_logZ"] is None:
        pars["gas_logZ"] = pars["logZ"]
    abunds = getNebAbunds(pars["set_name"],
                          pars["gas_logZ"],
                          dust=pars["dust"],
                          re_z=pars["re_z"])
    this_print('////////////////////////////////////')
    this_print('title {0}'.format(model_name.split('/')[-1]))
    this_print('////////////////////////////////////')
    this_print('set punch prefix "{0}"'.format(model_name))
    this_print('print line precision 6')
    ####
    if pars['par1'] == "age":
        pars['par1val'] = pars['age']
    if pars['par2'] == "logz":
        pars['par2val'] = pars['logZ']
        if pars['maxStellar'] is not None:
            if pars['logZ'] > pars['maxStellar']:
                pars['par2val'] = pars['maxStellar']
            else:
                pars['par2val'] = pars['logZ']
    this_print('table star "{0}" {1}={2:.2e} {3}={4:.2e}'.format(pars['cloudy_mod'], pars['par1'], pars['par1val'],pars['par2'], pars['par2val']))

    # logic: to add metallicity-dependent imbh at chosen fractional contribution (one frac contrib per run)
    if pars['table_sed']:
        table_sed_dict = {'0.01': 'M-2_00-qsoagn-i_45.sed', '0.05': 'M-2_19-qsoagn-i_45.sed', '0.1': 'M-3_28-qsoagn-i_45.sed', '0.2': 'M-3_96-qsoagn-i_45.sed', '0.4':  'M-4_59-qsoagn-i_45.sed', '0.7': 'M-5_24-qsoagn-i_45.sed', '1.0': 'M-5_89-qsoagn-i_45.sed'} # light seed
        #table_sed_dict = {'0.01': 'M-5_00-qsoagn-i_45.sed', '0.05': 'M-5_00-qsoagn-i_45.sed', '0.1': 'M-5_00-qsoagn-i_45.sed', '0.2': 'M-5_00-qsoagn-i_45.sed', '0.4':  'M-5_00-qsoagn-i_45.sed', '0.7': 'M-5_00-qsoagn-i_45.sed', '1.0': 'M-5_00-qsoagn-i_45.sed'} # if doing imbhFixedM05

        zsol_key = str(round(10**pars['logZ'],2))
        table_sed_file = table_sed_dict[zsol_key]

        if pars['table_sed_frac'] == 0.00:
            QH_tablesed = -30.0
            QH_stsed = pars['logQ']
        elif  pars['table_sed_frac'] == 1.00:
            QH_tablesed = pars['logQ']
            QH_stsed = -30.0
        else:
            QH_tablesed = np.log10(pars['table_sed_frac']*10**pars['logQ'])
            QH_stsed = np.log10((1-pars['table_sed_frac'])*10**pars['logQ'])
        this_print('Q(H) = {0:.3f} log'.format(QH_stsed))
        this_print('table sed "{0}"'.format(table_sed_file))
        this_print('Q(H) = {0:.3f} log'.format(QH_tablesed))

    elif pars['table_sed'] is None:
        pass
    if pars['use_Q'] and not pars['table_sed']:
        this_print('Q(H) = {0:.3f} log'.format(pars['logQ']))
    elif not pars['use_Q'] and not pars['table_sed']:
        this_print('ionization parameter = {0:.3f} log'.format(pars['logU']))


    ####
    this_print(abunds.solarstr)
    if pars['dust'] and pars['set_name'] != 'nicholls':
        #this_print('metals grains {0:.2f} log'.format(pars['gas_logZ']))
        this_print('grains {0:.2f} log'.format(pars['gas_logZ']))
    elif pars['dust'] and pars['set_name'] == 'nicholls':
        pass
    #else:
    #    #this_print('metals {0:.2f} log'.format(pars['gas_logZ']))
    for line in abunds.elem_strs:
        this_print(line)
    ####
    if pars['r_in_pc']:
        pc_to_cm = 3.08568e18
        r_out = np.log10(pars['r_inner']*pc_to_cm)
    else:
        r_out = pars['r_inner']
    if pars['use_extended_lines']:
        pkgdir = os.path.dirname(sys.modules['cloudyfsps'].__file__)
        linefile = 'cloudyLinesEXT.dat'
        #linefile =  os.path.join(pkgdir,'data/cloudyLinesEXT.dat')
        #linefile = pkg_resources.resource_filename(__name__,'data/cloudyLinesEXT.dat')
    else:
        pkgdir = os.path.dirname(sys.modules['cloudyfsps'].__file__)
        linefile = 'cloudyLines.dat'
        #linefile =  os.path.join(pkgdir,'data/cloudyLines.dat')
        #linefile = pkg_resources.resource_filename(__name__, 'data/cloudyLines.dat')
    this_print('radius {0:.3f} log'.format(r_out))
    this_print('hden {0:.3f} log'.format(np.log10(pars['dens'])))
    this_print('{}'.format(pars['geometry']))
    this_print('cosmic ray background')
    this_print('iterate to convergence max=5')
    #this_print('stop temperature 100.0') #changing stopping criteria to get partially ionized zone
    this_print('constant pressure') #changing stopping criteria to get partially ionized zone
    this_print("turbulence 2 km/s") # from chris's sims
    this_print("magnetic field -4") # from chris's sims
    this_print('stop efrac {0:.2f}'.format(pars['efrac']))
    this_print('save last linelist ".lin" "{}" absolute column'.format(linefile))
    this_print('save last outward continuum ".outwcont" units Angstrom no title')
    this_print('save last incident continuum ".inicont" units Angstrom no title')
    if len(pars["extras"]) > 0:
        this_print(pars["extras"])
    if pars["extra_output"]:
        this_print(extra_str)
    if pars["verbose"]:
        print("Input written in {0}".format(file_name))
        f.close()

extra_str = '''
save last radius ".rad"
save last physical conditions ".phys"
save last element hydrogen ".ele_H"
save last element helium ".ele_He"
save last element carbon ".ele_C"
save last element nitrogen ".ele_N"
save last element oxygen ".ele_O"
save last element sulphur ".ele_S"
save last element silicon ".ele_Si"
save last element iron ".ele_Fe"
save last hydrogen Lya ".H_lya"
save last hydrogen ionization ".H_ion"
save last lines emissivity ".emis"
H  1 6562.81A
H  1 4861.33A
H  1 4340.46A
H  1 4101.73A
He 1 3888.63A
He 1 4471.49A
He 1 5875.64A
He 2 1640.43A
He 2 4685.64A
O  1 6300.30A
O  1 6363.78A
O  2 3728.81A
O  2 3726.03A
O  3 5006.84A
Blnd 4363.00A
O  3 4958.91A
O  3 51.8004m
N  2 6583.45A
N  2 6548.05A
S  2 6730.82A
S  2 6716.44A
S  3 9068.62A
S  3 9530.62A
S  3 18.7078m
S  4 10.5076m
Ar 3 7135.79A
Ar 3 8.98898m
Ne 3 3868.76A
Ne 2 12.8101m
Ne 3 15.5509m
C  3 1908.73A
C  3 1906.68A
end of lines
'''

def writeMake(dir_=None):
    '''
    writes makefile that runs Cloudy on all files in directory with
    the same prefix.
    '''
    makefile = open("{0}/Makefile".format(dir_), "w")
    txt_exe = "CLOUDY = {0}\n".format(os.environ["CLOUDY_EXE"])
    txt = """
SRC = $(wildcard ${name}*.in)
OBJ = $(SRC:.in=.out)

#Usage: make -j n_proc name='NAME'
#optional: NAME is a generic name, all models named NAME*.in will be run

all: $(OBJ)

%.out: %.in
\t-$(CLOUDY) < $< > $@
#notice the previous line has TAB in first column
"""
    makefile.write(txt_exe)
    makefile.write(txt)
    makefile.close()

def runMake(dir_=None, n_proc=1, model_name=None):
    if dir_ is None:
        dir_="./"
    to_run = "cd {0} ; make -j {1:d}".format(dir_, n_proc)
    if model_name is not None:
        to_run += " name='{0}'".format(model_name)
    stdin = None
    stdout = subprocess.PIPE
    print("running: {0}".format(to_run))
    proc = subprocess.Popen(to_run, shell=True, stdout=stdout, stdin=stdin)
    proc.communicate()

def printParFile(dir_, mod_prefix, pars):
    '''
    prints parameter file for easy parsing later
    modnum, Z, a, U, R, logQ, n, efrac
    '''
    outfile = "{}{}.pars".format(dir_, mod_prefix)
    f = open(outfile, "w")
    for i in range(len(pars)):
        par = pars[i]
        if len(par) > 7:
            pstr = "{0} {1:.2f} {2:.2e} {3:.2f} {4:.2f} {5:.2f} {6:.2f} {7:.2f} {8:.2e}\n".format(i+1, *par)
        else:
            pstr = "{0} {1:.2f} {2:.2e} {3:.2f} {4:.2f} {5:.2f} {6:.2f} {7:.2f}\n".format(i+1, *par)
        f.write(pstr)
    f.close()
    return

def writeParamFiles(**kwargs):
    '''
    for making grids of parameters.
    can pass arrays of ages, logZs, logUs, nHs.
    cloudy_input.param_files(extras='extra line to add to input')
    '''
    nom_dict = {"dir_":"./output/",
                "model_prefix":"ZAU",
                "ages":np.arange(1.0e6, 6.0e6, 1.0e6),
                "logZs":np.linspace(-2.0, 0.2, 5),
                "logQs":np.linspace(45,49, 5),
                "logUs":np.linspace(-3.0, -1.0, 5),
                "r_inners":np.array([19.]),
                "nhs":np.array([100.0]),
                "pc_to_cm":False,
                "run_cloudy":False,
                "set_name":"dopita",
                "use_Q":True,
                "dust":True,
                "re_z":False,
                "cloudy_mod":"FSPS_SPS.mod",
                "verbose":False,
                "efracs":np.array([-1.0]),
                "geometry":"sphere",
                "write_makefile":False,
                "extras":"",
                "extra_output":False,
                "table_sed":False,
                "table_sed_frac": 0.0}
    for key, val in list(kwargs.items()):
        nom_dict[key] = val
    pars = kwargs.get("pars", None)
    if pars is None:
        print("{} ages, {} logZs, {} logUs".format(len(nom_dict["ages"]),
                                                   len(nom_dict["logZs"]),
                                                   len(nom_dict["logUs"])))
        pars = [(Z, a, U, R, calcForLogQ(logU=U, Rinner=10.0**R, nh=n), n, efrac) for Z in nom_dict["logZs"] for a in nom_dict["ages"] for U in nom_dict["logUs"] for R in nom_dict["r_inners"] for n in nom_dict["nhs"] for efrac in nom_dict["efracs"]]
    # Z, a, U, R, Q, n, efrac
    print("{} models".format(len(pars)))
    full_model_names = ["{}{}".format(nom_dict["model_prefix"], n+1)
                        for n in range(len(pars))]
    printParFile(nom_dict["dir_"], nom_dict["model_prefix"], pars)
    #--------------------------------------------
    for par, name in zip(pars, full_model_names):
        cloudyInput(nom_dict["dir_"],
                    name,
                    logZ=par[0],
                    age=par[1],
                    logU=par[2],
                    r_inner=par[3],
                    logQ=par[4],
                    dens=par[5],
                    efrac=par[6],
                    set_name=nom_dict["set_name"],
                    use_Q=nom_dict["use_Q"],
                    dust=nom_dict["dust"],
                    re_z=nom_dict["re_z"],
                    cloudy_mod=nom_dict["cloudy_mod"],
                    verbose=nom_dict["verbose"],
                    geometry=nom_dict["geometry"],
                    extras=nom_dict["extras"],
                    extra_output=nom_dict["extra_output"],
                    table_sed = nom_dict["table_sed"],
                    table_sed_frac = nom_dict["table_sed_frac"])
    #--------------------------------------------
    if nom_dict["write_makefile"]:
        writeMake(dir_=nom_dict["dir_"])
    #--------------------------------------------
    if nom_dict["run_cloudy"]:
        runMake(dir_=nom_dict["dir_"], n_proc=4, model_name=nom_dict["model_prefix"])
