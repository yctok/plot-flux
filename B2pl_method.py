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



def readProf(fname, wdir='.'):
    """
    reads contents of text file into two lists, returns them as numpy arrays
    """

    fname = path.join(wdir, fname)
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

    return x, y


def loadg(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    # read first line for case string and grid size
    line = lines[0]
    words = line.split()

    nw = int(words[-2])
    nh = int(words[-1])
    psi = np.linspace(0, 1, nw)

    # read in scalar parameters
    #   note: word size of 16 characters each is assumed for all of the data to be read

    # line 1
    line = lines[1]
    rdim = float(line[0:16])
    zdim = float(line[16:32])
    rcentr = float(line[32:48])
    rleft = float(line[48:64])
    zmid = float(line[64:80])

    # line 2
    line = lines[2]
    rmaxis = float(line[0:16])
    zmaxis = float(line[16:32])
    simag = float(line[32:48])
    sibry = float(line[48:64])
    bcentr = float(line[64:80])

    # line 3
    line = lines[3]
    current = float(line[0:16])

    # read in profiles
    #   start by reading entire file into single list then split into individual profiles
    #   first block has 5 profiles of length nw and one array of length nh*nw

    temp = []
    count = 0
    lnum = 5
    terms = 5 * nw + nw * nh
    while count < terms:
        line = lines[lnum]
        numchar = len(line)
        nwords = numchar // 16
        count1 = 0
        while count1 < nwords:
            i1 = count1 * 16
            i2 = i1 + 16
            temp.append(float(line[i1:i2]))
            count1 += 1
            count += 1
        lnum += 1

    fpol = temp[0:nw]
    pres = temp[nw:2 * nw]
    ffprime = temp[2 * nw:3 * nw]
    pprime = temp[3 * nw:4 * nw]
    psirz_temp = temp[4 * nw:(4 + nh) * nw]
    qpsi = temp[(4 + nh) * nw:]

    # split psirz up into 2D matrix
    count = 0
    psirz = []
    while count < nh:
        ind1 = count * nw
        ind2 = ind1 + nw
        psirz.append(psirz_temp[ind1:ind2])
        count += 1

    # scalars for length of boundary and limiter arrays
    line = lines[lnum]
    words = line.split()
    nbbbs = int(words[0])
    limitr = int(words[1])

    # read boundary and limiter points into temp array

    temp = []
    count = 0
    terms = 2 * (nbbbs + limitr)
    lnum += 1
    while count < terms:
        line = lines[lnum]
        numchar = len(line)
        nwords = numchar // 16
        count1 = 0
        while count1 < nwords:
            i1 = count1 * 16
            i2 = i1 + 16
            temp.append(float(line[i1:i2]))
            count1 += 1
            count += 1
        lnum += 1
    bdry_temp = temp[0:(2 * nbbbs)]
    limit_temp = temp[(2 * nbbbs):]

    # split boundary into (R,Z) pairs
    count = 0
    rbdry = []
    zbdry = []
    while count < len(bdry_temp) - 1:
        rbdry.append(bdry_temp[count])
        zbdry.append(bdry_temp[count + 1])
        count += 2

    # split limiter into (R,Z) pairs
    count = 0
    rlim = []
    zlim = []
    while count < len(limit_temp) - 1:
        rlim.append(limit_temp[count])
        zlim.append(limit_temp[count + 1])
        count += 2

    g = dict(nw=nw, nh=nh, rdim=rdim, zdim=zdim, rcentr=rcentr, rleft=rleft, zmid=zmid,
             rmaxis=rmaxis, zmaxis=zmaxis, simag=simag, sibry=sibry, current=current,
             fpol=np.array(fpol),
             ffprime=np.array(ffprime), pprime=np.array(pprime), psirz=np.array(psirz),
             qpsi=np.array(qpsi), nbbbs=nbbbs, bcentr=bcentr,
             pres=np.array(pres), limitr=limitr, rbdry=np.array(rbdry),
             zbdry=np.array(zbdry), rlim=np.array(rlim), zlim=np.array(zlim))

    return g


