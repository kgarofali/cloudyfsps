#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import object
import os

import numpy as np
from .generalTools import sym_to_name
from scipy.interpolate import InterpolatedUnivariateSpline as InterpUS

def getNebAbunds(set_name, logZ, dust=True, re_z=False, **kwargs):
    '''
    neb_abund.get_abunds(set_name, logZ, dust=True, re_z=False)
    set_name must be 'dopita', 'newdopita', 'cl01' or 'yeh'
    '''
    allowed_names = ['dopita', 'newdopita', 'cl01', 'yeh',
                     'varyNO', 'gutkin', 'UVbyler', 'varyCO', 'LIMS', 'nicholls']
    if set_name in allowed_names:
        return eval('{}({}, dust={}, re_z={})'.format(set_name, logZ, dust, re_z))
    else:
        raise IOError(allowed_names)

class abundSet(object):
    def __init__(self, set_name, logZ):
        '''
        overarching class for abundance sets.
        abundSet('dopita', 0.0)
        '''
        self.logZ = logZ
        self.abund_0 = load_abund(set_name)
        self.depl = load_depl(set_name)
        self.calcSpecial()
        self.calcFinal()
        self.inputStrings()

    def calcSpecial(self):
        return
    def calcFinal(self):
        return
    def inputStrings(self):
        self.solarstr = 'abundances {} {}'.format(self.solar, self.grains)
        elem_strs = []
        names = sym_to_name()
        for key in list(self.abund_0.keys()):
            elm = names[key]
            abund = self.__getattribute__(key)
            #if hasattr(self, 're_z'):
            #    if key != 'He':
            #        abund -= self.re_z
            outstr = 'element abundance {0} {1:.2f} log'.format(elm, abund)
            elem_strs.append(outstr)
        self.__setattr__('elem_strs', elem_strs)
        return

class dopita(abundSet):
    solar = 'old solar 84'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Dopita+2001: old solar abundances = 0.019
        ISM grains
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        if re_z:
            self.re_z = logZ
        else:
            self.re_z = 0.0
        abundSet.__init__(self, 'dopita', logZ)

    def calcSpecial(self):
        '''
        piece-wise function for nitrogen abund (step-function)
        functional form for helium
        '''
        def calc_N(logZ):
            if logZ <= -0.63:
                return -4.57 + logZ
            else:
                return -3.94 + (2.0*logZ)
        def calc_He(logZ):
            return np.log10(0.08096 + (0.02618*(10.0**logZ)))

        self.__setattr__('He', calc_He(self.logZ))
        self.__setattr__('N', calc_N(self.logZ)+self.depl['N'])
        return
    def calcFinal(self):
        '''
        apply depletions and scale with logZ
        '''
        [self.__setattr__(key, val+self.logZ+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return
class newdopita(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Abundances from Dopita (2013)
            Solar Abundances from Grevasse 2010 - z= 0.013
            includes smooth polynomial for N/O, C/O relationship
            functional form for He(z)
            new depletion factors
            ISM grains
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'newdopita', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024*(10.0**logZ)))
        def calc_CNO(logZ):
            oxy = np.array([7.39, 7.50, 7.69, 7.99, 8.17,
                    8.39, 8.69, 8.80, 8.99, 9.17, 9.39])
            nit = np.array([-6.61, -6.47, -6.23, -5.79, -5.51,
                    -5.14, -4.60, -4.40, -4.04, -3.67, -3.17])
            car = np.array([-5.58, -5.44, -5.20, -4.76, -4.48,
                    -4.11, -3.57, -3.37, -3.01, -2.64, -2.14])
            O = self.abund_0['O'] + logZ
            C = float(InterpUS(oxy, car, k=1)(O + 12.0))
            N = float(InterpUS(oxy, nit, k=1)(O + 12.0))
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val + self.depl[key])
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return
    def calcFinal(self):
        [self.__setattr__(key, val+self.logZ+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return
class UVbyler(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Abundances from Dopita (2013)
            Solar Abundances from Grevasse 2010 - z= 0.013
            New fit for N/O, C/O relationship
            functional form for He(z)
            new depletion factors
            ISM grains
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'UVbyler', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024*(10.0**logZ)))
        def calc_CNO(logZ):
            O = self.abund_0['O'] + logZ
            logOH = O + 12.0
            logCO = -0.8 + 0.14*(logOH - 8.0) + (0.192*np.log(1. + np.exp((logOH - 8.0)/0.2)))
            logNO = -1.5 + (0.1*np.log(1. + np.exp((logOH - 8.3)/0.1)))
            #C = np.log10((10.**O)*(10.**-0.789 + 10.**(4.105 + 1.263*O)))
            #N = np.log10((10.**O)*(10.**-1.579 + 10.**(3.579 + 1.526*O)))
            C = logCO + O
            N = logNO + O
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val + self.depl[key])
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return
    def calcFinal(self):
        [self.__setattr__(key, val+self.logZ+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return

class LIMS(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Solar UVByler abundances modified:
             Enhance alpha abundances +0.2 dex (O, Ne, Mg, Si, S, Ar, Ca)
                  Zhu+2010, Conroy+2014, Choi+2014
             Enhance C, N following PNe abundances (logNO ~ -0.5, logCO ~ 0)
                  Henry+2018, but also Karakas 2010, Maciel 2017
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'LIMS', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024*(10.0**logZ)))
        self.__setattr__('He', calc_He(self.logZ))
        return
    def calcFinal(self):
        [self.__setattr__(key, val+self.logZ+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return

class IIZw(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Abundances from Dopita (2013)
            Solar Abundances from Grevasse 2010 - z= 0.013
            New fit for N/O, C/O relationship
            functional form for He(z)
            new depletion factors
            ISM grains
            log O/H + 12. = 8.09
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'UVbyler', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024*(10.0**-0.6)))
        def calc_CNO(logZ):
            #O = self.abund_0['O'] + logZ
            #logOH = O + 12.0
            O = self.abund_0['O'] + -0.6
            logOH = 8.09
            #
            logCO = -0.8 + 0.14*(logOH - 8.0) +\
                    (0.192*np.log(1.+np.exp((logOH - 8.0)/0.2)))
            logNO = -1.5 + (0.1*np.log(1. + np.exp((logOH - 8.3)/0.1)))
            C = logCO + O + logZ
            N = logNO + O
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val)
         for key, val in zip(['C', 'O'], [C, O])]
        self.N = N + self.depl['N']
        return
    def calcFinal(self):
        [self.__setattr__(key, val+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return


class varyCO(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=0.0):
        '''
        arbitrarily vary C/O at fixed O.
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'UVbyler', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return np.log10(0.0737 + (0.024*(10.0**logZ)))
        def calc_CNO(logZ):
            O = self.abund_0['O'] + logZ
            logOH = O + 12.0
            logCO = -0.8 + 0.14*(logOH - 8.0) + (0.192*np.log(1. + np.exp((logOH - 8.0)/0.2)))
            logNO = -1.5 + (0.1*np.log(1. + np.exp((logOH - 8.3)/0.1)))
            #C = np.log10((10.**O)*(10.**-0.789 + 10.**(4.105 + 1.263*O)))
            #N = np.log10((10.**O)*(10.**-1.579 + 10.**(3.579 + 1.526*O)))
            C = logCO + O
            N = logNO + O
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        self.__setattr__('C', C + self.depl['C'] + self.re_z)
        self.__setattr__('N', N + self.depl['N'])
        self.__setattr__('O', O + self.depl['O'])
        return
    def calcFinal(self):
        [self.__setattr__(key, val + self.re_z + self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return

class gutkin(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        Gutkin+2016
            PARSEC metallicity (Bressan+2012)
            based on Grevesse+Sauvel (1998) and Caffau+2011
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'gutkin', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            Z = (10.**logZ)*0.01524
            Y = 0.2485 + 1.7756*Z
            X = 1. - Y - Z
            return np.log10(Y/X/4.)
        def calc_CNO(logZ):
            O = self.abund_0['O'] + logZ
            N = np.log10((0.41 * 10.**O)*(10.**-1.6 + 10.**(2.33 + O)))
            C = self.abund_0['C'] + logZ
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val)
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return
    def calcFinal(self):
        [self.__setattr__(key, val)
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return

class varyNO(abundSet):
    solar = 'GASS10'
    def __init__(self, logZ, dust=True, re_z=False):
        '''
        varying N at fixed O.
        '''
        if dust:
            self.grains = 'no grains\ngrains ISM'
        else:
            self.grains = 'no grains'
        self.re_z=re_z
        abundSet.__init__(self, 'dopita', logZ)

    def calcSpecial(self):
        def calc_He(logZ):
            return -1.01
        def calc_CNO(logZ):
            oxy = np.array([7.39, 7.50, 7.69, 7.99, 8.17,
                    8.39, 8.69, 8.80, 8.99, 9.17, 9.39])
            nit = np.array([-6.61, -6.47, -6.23, -5.79, -5.51,
                    -5.14, -4.60, -4.40, -4.04, -3.67, -3.17])
            car = np.array([-5.58, -5.44, -5.20, -4.76, -4.48,
                    -4.11, -3.57, -3.37, -3.01, -2.64, -2.14])
            O = self.abund_0['O']
            C = float(InterpUS(oxy, car, k=1)(O + 12.0))
            N = float(InterpUS(oxy, nit, k=1)(O + logZ + 12.0))
            return C, N, O
        self.__setattr__('He', calc_He(self.logZ))
        C, N, O = calc_CNO(self.logZ)
        [self.__setattr__(key, val)
         for key, val in zip(['C', 'N', 'O'], [C, N, O])]
        return
    def calcFinal(self):
        [self.__setattr__(key, val)
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return

class nicholls(abundSet):

    def __init__(self, logZ, dust=True, re_z=False):
        abund_logZ = {'0.0002': 'GC_ZO_0010.abn', '0.001': 'GC_ZO_0050.abn', '0.002': 'GC_ZO_0100.abn', '0.003': 'GC_ZO_0150.abn', '0.004': 'GC_ZO_0200.abn', '0.006': 'GC_ZO_0300.abn', '0.008': 'GC_ZO_0400.abn', '0.01': 'GC_ZO_0500.abn', '0.014': 'GC_ZO_0700.abn', '0.02': 'GC_ZO_1000.abn', '0.03': 'GC_ZO_1500.abn', '0.04': 'GC_ZO_2000.abn'}
        abund_key = str(round(round(10**logZ,4)*0.02,4))
        abund_file = os.environ['CLOUDY_DATA_PATH']+'abundances/'+str(abund_logZ[abund_key])
        if dust:
            D_G_logZ = {'0.0002': 0.00001217, '0.001': 0.00179, '0.002': 0.0154, '0.003': 0.0542, '0.004': 0.132, '0.006': 0.332, '0.008': 0.442, '0.01': 0.555, '0.014': 0.785, '0.02': 1.136, '0.03': 1.74, '0.04': 2.382}
            D_G_key = str(round(round(10**logZ,4)*0.02,4))
            ograins = '\ngrains orion '+str(D_G_logZ[D_G_key])
            pahgrains = '\ngrains PAH '+str(D_G_logZ[D_G_key])
            self.grains = ograins+pahgrains

        else:
            self.grains = 'no grains'
        self.re_z=re_z
        self.solar = '"'+str(abund_logZ[abund_key])+'"'
        abundSet.__init__(self, 'nicholls', logZ)

    def calcSpecial(self):
        pass
    def calcFinal(self):

        abund_logZ = {'0.0002': 'GC_ZO_0010.abn', '0.001': 'GC_ZO_0050.abn', '0.002': 'GC_ZO_0100.abn', '0.003': 'GC_ZO_0150.abn', '0.004': 'GC_ZO_0200.abn', '0.006': 'GC_ZO_0300.abn', '0.008': 'GC_ZO_0400.abn', '0.01': 'GC_ZO_0500.abn', '0.014': 'GC_ZO_0700.abn', '0.02': 'GC_ZO_1000.abn', '0.03': 'GC_ZO_1500.abn', '0.04': 'GC_ZO_2000.abn'}
        elemd = {'Helium': 'He', 'Lithium': 'Li', 'Beryllium': 'Be', 'Boron': 'B', 'Carbon': 'C', 'Nitrogen': 'N', 'Oxygen': 'O', 'Fluorine': 'F', 'Neon': 'Ne', 'Sodium': 'Na', 'Magnesium': 'Mg', 'Aluminium': 'Al', 'Silicon': 'Si', 'Phosphorus': 'P', 'Sulphur': 'S', 'Chlorine': 'Cl', 'Argon': 'Ar', 'Potassium': 'K', 'Calcium': 'Ca', 'Scandium': 'Sc', 'Titanium': 'Ti', 'Vanadium': 'V', 'Chromium': 'Cr', 'Manganese': 'Mn', 'Iron': 'Fe', 'Cobalt': 'Co', 'Nickel': 'Ni', 'Copper': 'Cu', 'Zinc': 'Zn'}
        abund_key = str(round(round(10**self.logZ,4)*0.02,4))
        abund_file = os.environ['CLOUDY_DATA_PATH']+'abundances/'+str(abund_logZ[abund_key])
        with open(abund_file) as f:
            for line in f:
                if not line.startswith('#'):
                    (key, val) = line.split()
                    if key == 'Hydrogen':
                        pass
                    else:
                        self.abund_0[elemd[key]] = np.log10(float(val))
        # nell's would take dictionary below (at solar) and add logZ value to scale with metallicity. maybe this is equivalent to using chris's different abund files as above. here we just use the actual file values, then add depletion factors, and don't add logZ (i.e, no "val+self.logZ")
        [self.__setattr__(key, val+self.depl[key])
         for key, val in list(self.abund_0.items()) if not hasattr(self, key)]
        return


def load_abund(set_name):
    if set_name == 'dopita':
        adict = dict(He=-1.01,
                     C=-3.44,
                     N=-3.95,
                     O=-3.07,
                     Ne=-3.91,
                     Mg=-4.42,
                     Si=-4.45,
                     S=-4.79,
                     Ar=-5.44,
                     Ca=-5.64,
                     Fe=-4.33,
                     F=-7.52,
                     Na=-5.69,
                     Al=-5.53,
                     P=-6.43,
                     Cl=-6.73,
                     K=-6.87,
                     Ti=-6.96,
                     Cr=-6.32,
                     Mn=-6.47,
                     Co=-7.08,
                     Ni=-5.75,
                     Cu=-7.73,
                     Zn=-7.34)
    elif set_name == 'newdopita':
        adict = dict(He=-1.01,
                     C=-3.57,
                     N=-4.60,
                     O=-3.31,
                     Ne=-4.07,
                     Na=-5.75,
                     Mg=-4.40,
                     Al=-5.55,
                     Si=-4.49,
                     S=-4.86,
                     Cl=-6.63,
                     Ar=-5.60,
                     Ca=-5.66,
                     Fe=-4.50,
                     Ni=-5.78,
                     F=-7.44,
                     P=-6.59,
                     K=-6.97,
                     Cr=-6.36,
                     Ti=-7.05,
                     Mn=-6.57,
                     Co=-7.01,
                     Cu=-7.81,
                     Zn=-7.44)
    elif set_name == 'UVbyler':
        adict = dict(He=-1.01,
                     C=-3.57,
                     N=-4.17,
                     O=-3.31,
                     Ne=-4.07,
                     Na=-5.75,
                     Mg=-4.40,
                     Al=-5.55,
                     Si=-4.49,
                     S=-4.86,
                     Cl=-6.63,
                     Ar=-5.60,
                     Ca=-5.66,
                     Fe=-4.50,
                     Ni=-5.78,
                     F=-7.44,
                     P=-6.59,
                     K=-6.97,
                     Cr=-6.36,
                     Ti=-7.05,
                     Mn=-6.57,
                     Co=-7.01,
                     Cu=-7.81,
                     Zn=-7.44)
    elif set_name == 'LIMS':
        adict = dict(He=-1.01,
                     C=-3.11,
                     N=-3.61,
                     O=-3.11,
                     Ne=-3.87,
                     Na=-5.75,
                     Mg=-4.20,
                     Al=-5.55,
                     Si=-4.29,
                     S=-4.66,
                     Cl=-6.63,
                     Ar=-5.40,
                     Ca=-5.46,
                     Fe=-4.50,
                     Ni=-5.78,
                     F=-7.44,
                     P=-6.59,
                     K=-6.97,
                     Cr=-6.36,
                     Ti=-7.05,
                     Mn=-6.57,
                     Co=-7.01,
                     Cu=-7.81,
                     Zn=-7.44)
    elif set_name == 'gutkin':
        adict = dict(He=-1.01,
                     C=-3.53,
                     N=-4.32,
                     O=-3.17,
                     F=-7.47,
                     Ne=-4.01,
                     Na=-5.70,
                     Mg=-4.45,
                     Al=-5.56,
                     Si=-4.48,
                     P=-6.57,
                     S=-4.87,
                     Cl=-6.53,
                     Ar=-5.63,
                     K=-6.92,
                     Ca=-5.67,
                     Sc=-8.86,
                     Ti=-7.01,
                     V=-8.03,
                     Cr=-6.36,
                     Mn=-6.64,
                     Fe=-4.51,
                     Co=-7.11,
                     Ni=-5.78,
                     Cu=-7.82,
                     Zn=-7.43)
    elif set_name =='nicholls':
        # from table 1 in richardson+2022
        adict = dict(He=-1.01,
                     Li=-8.722,
                     Be=-10.68,
                     B=-9.193,
                     C=-3.577,
                     N=-4.21,
                     O=-3.24,
                     F=-7.56,
                     Ne=-3.91,
                     Na=-5.79,
                     Mg=-4.44,
                     Al=-5.57,
                     Si=-4.50,
                     P=-6.59,
                     S=-4.88,
                     Cl=-6.75,
                     Ar=-5.60,
                     K=-6.96,
                     Ca=-5.68,
                     Sc=-8.84,
                     Ti=-7.07,
                     V=-8.11,
                     Cr=-6.38,
                     Mn=-6.58,
                     Fe=-4.48,
                     Co=-7.07,
                     Ni=-5.80,
                     Cu=-7.82,
                     Zn=-7.44)
    return adict

def load_depl(set_name):
    if set_name == 'dopita':
        ddict = dict(C=-0.30,
                     N=-0.22,
                     O=-0.22,
                     Ne=0.0,
                     Mg=-0.70,
                     Si=-1.0,
                     S=0.0,
                     Ar=0.0,
                     Ca=-2.52,
                     Fe=-2.0,
                     F=0.0,
                     Na=0.0,
                     Al=0.0,
                     P=0.0,
                     Cl=0.0,
                     K=0.0,
                     Ti=0.0,
                     Cr=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Ni=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'newdopita':
        ddict = dict(He=0.00,
                     C=-0.30,
                     N=-0.05,
                     O=-0.07,
                     Ne=0.00,
                     Na=-1.00,
                     Mg=-1.08,
                     Al=-1.39,
                     Si=-0.81,
                     S=0.00,
                     Cl=-1.00,
                     Ar=0.00,
                     Ca=-2.52,
                     Fe=-1.31,
                     Ni=-2.00,
                     F=0.0,
                     P=0.0,
                     K=0.0,
                     Cr=0.0,
                     Ti=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'UVbyler':
        ddict = dict(He=0.00,
                     C=-0.30,
                     N=-0.05,
                     O=-0.07,
                     Ne=0.00,
                     Na=-1.00,
                     Mg=-1.08,
                     Al=-1.39,
                     Si=-0.81,
                     S=0.00,
                     Cl=-1.00,
                     Ar=0.00,
                     Ca=-2.52,
                     Fe=-1.31,
                     Ni=-2.00,
                     F=0.0,
                     P=0.0,
                     K=0.0,
                     Cr=0.0,
                     Ti=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'LIMS':
        ddict = dict(He=0.00,
                     C=-0.30,
                     N=-0.05,
                     O=-0.07,
                     Ne=0.00,
                     Na=-1.00,
                     Mg=-1.08,
                     Al=-1.39,
                     Si=-0.81,
                     S=0.00,
                     Cl=-1.00,
                     Ar=0.00,
                     Ca=-2.52,
                     Fe=-1.31,
                     Ni=-2.00,
                     F=0.0,
                     P=0.0,
                     K=0.0,
                     Cr=0.0,
                     Ti=0.0,
                     Mn=0.0,
                     Co=0.0,
                     Cu=0.0,
                     Zn=0.0)
    elif set_name == 'gutkin':
        ddict = dict(He=0.00,
                     Li=-0.8,
                     C=-0.30,
                     O=-0.15,
                     Na=-0.60,
                     Mg=-0.70,
                     Al=-1.70,
                     Si=-1.00,
                     Cl=-0.30,
                     Ca=-2.52,
                     Fe=-2.00,
                     Ni=-1.40)
    elif set_name =='nicholls':
        # from table 1 in richardson+2022
        ddict = dict(He=0.00,
                     Li=-0.524,
                     Be=-0.274,
                     B=-0.546,
                     C=-0.120,
                     N=0.00,
                     O=-0.112,
                     F=-0.147,
                     Ne=0.00,
                     Na=-0.638,
                     Mg=-0.659,
                     Al=-1.602,
                     Si=-0.625,
                     P=0.00,
                     S=0.00,
                     Cl=-0.037,
                     Ar=0.00,
                     K=-0.614,
                     Ca=-2.398,
                     Sc=-1.533,
                     Ti=-1.928,
                     V=-1.159,
                     Cr=-1.379,
                     Mn=-1.134,
                     Fe=-1.510,
                     Co=-1.343,
                     Ni=-1.517,
                     Cu=-0.757,
                     Zn=-0.075)
    return ddict
