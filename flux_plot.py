# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 15:45:19 2024

@author: ychuang
"""


import os
import numpy as np
import matplotlib.pyplot as plt

import B2pl_method as bm

plt.rcParams.update({'font.weight': 'bold'})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.facecolor':'w'})
plt.rcParams.update({'mathtext.default': 'regular'})

eV = 1.60217662e-19


class fluxplot:

    def __init__(self, workdir, gfile_loc, b2plot_dev_x11 = False):
        """
        Make sure you source setup.csh or setup.ksh before running any of this!
        
        Inputs:
          workdir         Directory with the SOLPS outputs
          gfile_loc       location of corresponding g file
          impurity_list   List of all the impurity species included in the plasma simulation
          b2plot_dev_x11  Set to True if you want a figure to pop up for every b2plot call
        """
        shot_loc_in_gfile_string = gfile_loc.rfind('g')
        shot = int(gfile_loc[shot_loc_in_gfile_string+1 : shot_loc_in_gfile_string+7])

        shot_ind_in_workdir = workdir.rfind(str(shot))
        if shot_ind_in_workdir > 0:
            workdir_short = workdir[shot_ind_in_workdir:]
        else:
            workdir_short = None

        self.data = {'workdir':workdir, 'workdir_short':workdir_short, 'gfile_loc': gfile_loc,
                     'expData':{'fitProfs':{}}, 'solpsData':{'profiles':{}}}

        if 'B2PLOT_DEV' in os.environ.keys():
            if (not b2plot_dev_x11) and os.environ['B2PLOT_DEV'] == 'x11 ps':
                print("Changing environment variable B2PLOT_DEV to 'ps'")
                os.environ['B2PLOT_DEV'] = 'ps'
        else:
            print('WARNING: Need to source setup.csh for SOLPS-ITER distribution for complete SOLPSxport workflow')

        self.timeid = None

    

    def calcPsiVals(self, plotit = False, shift=1):
        """
        Call b2plot to get the locations of each grid cell in psin space

        Saves the values to dictionaries in self.data['solpsData']
        """
        from scipy import interpolate

        """
        Find grid corners first:
          0: lower left
          1: lower right
          2: upper left
          3: upper right

        Average location of cells 0 and 2 for middle of 'top' surface, 
        which is the top looking at outboard midplane
        Don't average over whole cell, dR << dZ at outboard midplane 
        and surface has curvature, so psin will be low

        jxa = poloidal cell index for the outer midplane
        crx = radial coordinate corner of grid [m]
        cry = vertical coordinate corner of grid [m]
        writ = write b2plot.write file
        f.y = plot against y
        """
        wdir = self.data['workdir']

        dsa, crLowerLeft = bm.B2pl('0 crx writ jxa f.y', wdir = wdir)
        # dummy, crLowerRight = B2pl('1 crx writ jxa f.y', wdir = wdir)
        # Only 2 unique psi values per cell, grab 0 and 2
        dummy, crUpperLeft = bm.B2pl('2 crx writ jxa f.y', wdir = wdir)  # all x inds are the same
        dummy, czLowerLeft = bm.B2pl('0 cry writ jxa f.y', wdir = wdir)
        dummy, czUpperLeft = bm.B2pl('2 cry writ jxa f.y', wdir = wdir)
        ncells = len(dummy)

        g = bm.loadg(self.data['gfile_loc'])
        d = float(shift)
        psiN = (g['psirz'] - g['simag']) / (g['sibry'] - g['simag'])

        dR = g['rdim'] / (g['nw'] - 1)
        dZ = g['zdim'] / (g['nh'] - 1)

        gR = []
        for i in range(g['nw']):
            gR.append(g['rleft'] + i * dR + d)

        gZ = []
        for i in range(g['nh']):
            gZ.append(g['zmid'] - 0.5 * g['zdim'] + i * dZ)

        gR = np.array(gR)
        gZ = np.array(gZ)

        R_solps_top = 0.5 * (np.array(crLowerLeft) + np.array(crUpperLeft))
        Z_solps_top = 0.5 * (np.array(czLowerLeft) + np.array(czUpperLeft))

        psiNinterp = interpolate.interp2d(gR, gZ, psiN, kind = 'cubic')

        psi_solps = np.zeros(ncells)
        for i in range(ncells):
            psi_solps_LL = psiNinterp(crLowerLeft[i], czLowerLeft[i])
            psi_solps_UL = psiNinterp(crUpperLeft[i], czUpperLeft[i])
            psi_solps[i] = np.mean([psi_solps_LL,psi_solps_UL])
        
        

        self.data['solpsData']['crLowerLeft'] = np.array(crLowerLeft)
        self.data['solpsData']['czLowerLeft'] = np.array(czLowerLeft)
        self.data['solpsData']['dsa'] = np.array(dsa)
        self.data['solpsData']['psiSOLPS'] = np.array(psi_solps)
        
        #from IPython import embed; embed()

        if plotit:
            psiN_range = [np.min(psi_solps), np.max(psi_solps)]

            psiN_copy = psiN.copy()
            psiN_copy[psiN > psiN_range[1]] = np.nan
            psiN_copy[psiN < psiN_range[0]] = np.nan
            psin_masked = np.ma.masked_invalid(psiN_copy)

            plt.figure()
            plt.contourf(gR, gZ, psin_masked, levels=[psiN_range[0], 1, psiN_range[1]], colors = [])
            # plt.pcolormesh(gR, gZ, psin_masked, cmap = 'inferno')
            # plt.colorbar(ticks = [0.25,0.5,0.75,1])
            plt.plot(g['rlim'], g['zlim'], 'k', lw = 2)
            # plt.contour(gR, gZ, psiN, [1], colors='k')
            # plt.contour(gR, gZ, psiN, [psiN_range[0]], colors='r', linestyles='dashed')
            # plt.contour(gR, gZ, psiN, [psiN_range[1]], colors='r', linestyles='dashed')
            gfile_name = self.data['gfile_loc'][self.data['gfile_loc'].rfind('/')+1:]
            plt.title(gfile_name)
            plt.axis('equal')
            plt.xlabel('R (m)')
            plt.ylabel('Z (m)')
            plt.xticks([1.0, 1.5, 2.0, 2.5])
            plt.yticks(np.arange(-1.5, 1.6, 0.5))
            plt.xlim([np.min(gR), np.max(gR)])
            plt.ylim([np.min(gZ), np.max(gZ)])
            plt.plot(R_solps_top, Z_solps_top, 'g', lw=3)
            plt.plot([1.94, 1.94], [-1.5, 1.5], ':k', lw=1)  # Thomson laser path
            plt.tight_layout()

            plt.figure()
            plt.plot(R_solps_top, psi_solps, 'k', lw = 2)
            plt.xlabel('R at midplane (m)')
            plt.ylabel('$\psi_N$')

            plt.show(block = False)






    def getSOLPSlast10Profs(self, plotit = False, use_existing_last10 = False):
        """
        Generates and reads the .last10 files (produced by running '2d_profiles', which looks
        at the last 10 time steps in the run.log file)
        """
        working = str(self.data['workdir'])

        olddir = os.getcwd()
        os.chdir(working)
        if working[-1] != '/': working += '/'
        # Call 2d_profiles as default, so you don't accidentally look at old time steps
        if (not use_existing_last10) or (not os.path.isfile(working + 'ne3da.last10')):
            print("Calling '2d_profiles' in directory: " + working)
            os.system('2d_profiles')

        rx, ne_ = bm.readProf('ne3da.last10')
        rx, dn_ = bm.readProf('dn3da.last10')
        rx, te_ = bm.readProf('te3da.last10')
        rx, ke_ = bm.readProf('ke3da.last10')
        rx, ti_ = bm.readProf('ti3da.last10')
        rx, ki_ = bm.readProf('ki3da.last10')
        
        os.chdir(olddir)

        # Cast everything as np array so it doesn't break later when performing math operations
        
        ne = np.array(ne_)
        dn = np.array(dn_)
        te = np.array(te_)
        ke = np.array(ke_)
        ti = np.array(ti_)
        ki = np.array(ki_)

        last10_dic = {'rx':rx,'ne':ne,'dn':dn,'te':te,'ke':ke,'ti':ti,'ki':ki}
    
        self.data['solpsData']['last10'] = last10_dic
        
        if plotit:
            if 'psiSOLPS' in self.data['solpsData'].keys():
                psi = self.data['solpsData']['psiSOLPS']

                f, ax = plt.subplots(2, 3, sharex = 'all')
                ax[0, 0].plot(psi, ne / 1.0e19, 'r', lw = 2, label = 'SOLPS')
                ax[0, 0].set_ylabel('n$_e$ (10$^{19}$ m$^{-3}$)')
                ax[0, 0].grid('on')

                ax[1, 0].plot(psi, dn, '-ok', lw = 2)
                ax[1, 0].set_ylabel('D')
                ax[1, 0].set_xlabel('$\psi_N$')

                ax[0, 1].plot(psi, te / 1.0e3, 'r', lw = 2, label = 'SOLPS')
                ax[0, 1].set_ylabel('Te (keV)')

                ax[1, 1].plot(psi, ke, 'b', lw = 2)
                ax[1, 1].set_ylabel('$\chi_e$')
                ax[1, 1].set_xlabel('$\psi_N$')
                ax[1, 1].set_xlim([np.min(psi) - 0.01, np.max(psi) + 0.01])

                ax[0, 2].plot(psi, ti / 1.0e3, 'r', lw = 2, label = 'SOLPS')
                ax[0, 2].set_ylabel('Ti (keV)')

                ax[1, 2].plot(psi, ki, 'b', lw = 2)
                ax[1, 2].set_ylabel('$\chi_i$')
                ax[1, 2].set_xlabel('$\psi_N$')
                ax[0, 0].set_title('last10 profiles')
                plt.tight_layout()
                
                for i in range(2):
                    for j in range(3):
                        ax[i,j].grid('on')
                        
                plt.figure()
                plt.plot(psi, dn / ke, 'k', lw=3, label = 'D / $\chi_e$')
                plt.plot(psi, dn / ki, 'r', lw=3, label = 'D / $\chi_i$')
                plt.grid('on')
                ax[0, 0].set_xticks(np.arange(0.84, 1.05, 0.04))
                plt.legend(loc='best')
                        
            else:
                f, ax = plt.subplots(3, sharex = 'all')
                ax[0].plot(rx, ne*1e-19, '-kx', lw=2)
                ax[0].set_ylabel('n$_e$ (10$^{19}$ m$^{-3}$)')
                ax[0].grid('on')
                
                ax[1].plot(rx, te, '-rx', lw=2, label = 'Te')
                ax[1].plot(rx, ti, '-bx', lw=2, label = 'Ti')
                ax[1].set_ylabel('T (eV)')
                ax[1].legend(loc='best')
                ax[1].grid('on')

                ax[2].plot(rx, dn, '-kx', lw=3, label = 'dn')
                ax[2].plot(rx, ke, '-gx', lw=3, label = 'ke')
                ax[2].plot(rx, ki, '-mx', lw=1, label = 'ki')
                ax[2].legend(loc='best')
                ax[2].set_ylabel('D or $\chi$')
                plt.grid('on')
                plt.tight_layout()
    
                ax[-1].set_xlabel('rx')
                
            plt.show(block = False)
    
      

    def getSOLPSfluxProfs(self, plotit):
        
        
        """
        Calls b2plot to get the particle flux profiles
        """
        # x variable is identical for all of these
        x_fTot, fluxTot = bm.B2pl("fnay za m* 0 0 sumz sy m/ writ jxa f.y")
        x_fTot, fluxD = bm.B2pl("fnay 1 zsel sy m/ writ jxa f.y")
        dummy, fluxConv = bm.B2pl("na za m* vlay m* 0 0 sumz writ jxa f.y")
        dummy, na = bm.B2pl("na 0 0 sumz writ jxa f.y")
        # dummy, hy1 = bm.B2pl("hy1 writ jxa f.y")  # not used anymore
        dummy, qe = bm.B2pl("fhey sy m/ writ jxa f.y")
        dummy, qi = bm.B2pl("fhiy sy m/ writ jxa f.y")
        dummy, core_qe = bm.B2pl("fhey sy m/ writ 0 f.x")
        dummy, core_qe_dat = bm.B2pl("fhey writ 0 f.x")
        dummy, sy = bm.B2pl("sy writ 0 f.x")
    
        for c in [fluxTot, fluxConv]:
            if not c:
                print("WARNING: Variable not populated by b2plot in getSOLPSfluxProfs")
                print("  Make sure ncl_ncar and netcdf modules are loaded")
                break
    
        self.data['solpsData']['profiles']['x_fTot'] = np.array(x_fTot)
        self.data['solpsData']['profiles']['fluxTot'] = np.array(fluxTot)
        self.data['solpsData']['profiles']['fluxD'] = np.array(fluxD)
        self.data['solpsData']['profiles']['fluxConv'] = np.array(fluxConv)
        self.data['solpsData']['profiles']['na'] = np.array(na)
        # self.data['solpsData']['profiles']['hy1'] = np.array(hy1)  # not used anymore
        self.data['solpsData']['profiles']['qe'] = np.array(qe)
        self.data['solpsData']['profiles']['qi'] = np.array(qi)
        
        print('the fluxD is:')
        print(fluxD)
        print('the fluxConv is:')
        print(fluxConv)
        print('the qe is:')
        print(qe)
        print('the radial cell area is:')
        print(sy)
        # print('the core_qe is:')
        # print(core_qe)
        # print('the core_qe_dat is:')
        # print(core_qe_dat)
    
        if plotit:
                
            # Check electron density from last10 profs for consistency
            ne_last10 = self.data['solpsData']['last10']['ne']
            rx_last10 = self.data['solpsData']['last10']['rx']  # very slightly different...
    
            f, ax = plt.subplots(3, sharex = 'all')
    
            ax[0].plot(rx_last10, ne_last10, '-kx', lw = 1, label = 'ne_last10')
            ax[0].plot(x_fTot, na, '--r*', lw=2, label = 'na')
            ax[0].set_ylabel('n (m$^{-3}$)')
            ax[0].legend(loc='best')
            ax[0].grid('on')
            
            """
            if self.data['workdir_short'] is not None:
                ax[0].set_title(self.data['workdir_short'], fontsize=10)
            else:
                ax[0].set_title('DIII-D shot ' + str(self.data['shot']) +
                                ', ' + str(self.timeid) + ' ms')
                
            """
    
            ax[1].plot(x_fTot, fluxTot, '-ko', lw = 2, label = 'Tot')
            ax[1].legend(loc='best')
            ax[1].set_ylabel('$\Gamma_tot$')
            ax[1].grid('on')
            
            

            ax[2].plot(x_fTot, fluxD, '-bx', lw = 2, label = 'Conv')
            ax[2].legend(loc='best')
            ax[2].set_ylabel('$\Gamma_D$')
            ax[2].grid('on')
            
            
            
            ax[-1].set_xlabel('x')
            
            ax[0].set_xlim([np.min(x_fTot) - 0.01, np.max(x_fTot) + 0.01])
            plt.show(block = True)