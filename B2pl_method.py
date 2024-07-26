# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 16:20:53 2024

@author: ychuang
"""




from os import path, environ
import numpy as np




def set_b2plot_dev(verbose = False):
    if 'B2PLOT_DEV' in environ.keys():
        if environ['B2PLOT_DEV'] == 'ps':
            if verbose:
                print('B2PLOT_DEV already set to ps')
        else:
            print("Changing environment variable B2PLOT_DEV to 'ps' for B2plot calls")
            environ['B2PLOT_DEV'] = 'ps'
    else:
        print('WARNING: Need to source setup.csh for a SOLPS-ITER distribution to enable B2plot calls')

        



def B2pl(cmds, wdir='.', debug=False):
    # import sys
    import subprocess
    """
    runs B2plot with the commands used in the call and reads contents of the resulting
    b2plot.write file into two lists
    
    ** Make sure you've sourced the setup script first, or this won't work! **
    **  Make sure B2PLOT_DEV is set to 'ps'
    """

    if debug:
        cmdstr = 'echo "' + cmds + '" | b2plot'
        print(cmdstr)
    else:
        # cmdstr = 'echo "' + cmds + '" | b2plot'
        # cmdstr = 'echo "' + cmds + '" | b2plot >&/dev/null'
        cmdstr = 'echo "' + cmds + '" | b2plot 2>/dev/null'
        testcmd = subprocess.check_output(cmdstr,shell=True)
        if testcmd == b'':
            print('\nERROR: b2plot command was not successful, is the case still running?')
            print('Command was: ',cmdstr)
            raise OSError
            
    fname = path.join(wdir, 'b2pl.exe.dir', 'b2plot.write')
    if not path.exists(fname):
        print('B2Plot writing failed for call:')
        print(cmds)
        print('in directory: ' + wdir + '\n')
        raise OSError

    # from IPython import embed; embed()

    x, y = [], []
    with open(fname) as f:
        lines = f.readlines()

    for line in lines:
        elements = line.split()
        if elements[0] == '#':
            pass
        else:
            x.append(float(elements[0]))
            y.append(float(elements[1]))
    x = x[0:(len(x) // 2)]  # used to be: x=x[0:(len(x)/2)-1], chopped final value
    y = y[0:(len(y) // 2)]

    return x, y        