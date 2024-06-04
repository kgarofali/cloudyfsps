#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
from __future__ import unicode_literals

import os
import sys
import subprocess
import numpy as np
from cloudyfsps.cloudyOutputTools import formatCloudyOutput

use_extended_lines=True
write_line_lum=True

def check_warning_ok(fl, mod_prefix, modnum):
    with open(fl) as of:
        if ('WARNING:' in of.read()):
            if ('OK' in of.read()):
                return mod_prefix+modnum+", F, warning, OK"
            else:
                 return mod_prefix+modnum+", F, warning, not OK"
        else:
            return mod_prefix+modnum+", T, OK"

def get_pars(dir_, mod_prefix, modnum):
    '''
    reads from your/mod/dir/PREFIX.pars
    order is:
    logZ, Age, logU, logR, logQ, nH, efrac, gas_logZ
    '''
    data = np.genfromtxt(dir_+mod_prefix+".pars")
    return data[np.int(modnum)-1, 1:]

def output_OK(fl, dir_, mod_prefix, modnum):
    '''
    looks in PREFIX123.out for 'Cloudy exited OK'
    '''
    ecode_file = dir_+mod_prefix+".ecode"
    with open(ecode_file, 'a') as ecode:
        ecode.write(check_warning_ok(fl, mod_prefix, modnum)+"\n")

    f = open(fl, 'r')
    ll = f.readlines()[-1]

    if 'OK' in ll:
        return True
    else:
        return False
def main(argv):
    fulloutfile=argv[0]
    dir_='/'.join(fulloutfile.split('/')[0:-1])+'/'
    fname = fulloutfile.split('/')[-1]
    mod_prefix=fname.split('.out')[0][0:3]
    modnum=fname.split('.out')[0][3::]
    if os.path.exists(fulloutfile):
        if output_OK(fulloutfile, dir_, mod_prefix, modnum):
            print("Cloudy exited OK.")
            print("Formatting output...")
            modpars = get_pars(dir_, mod_prefix, modnum)
            formatCloudyOutput(dir_, mod_prefix, modnum, modpars,
                               use_extended_lines=use_extended_lines)
        else:
            print("CLOUDY ERROR in {}. Stopping.".format(fname))
    else:
        print("ERROR, no outfile. Stopping.")

if __name__ == "__main__":
    main(sys.argv[1:])
