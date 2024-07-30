# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 16:44:40 2024

@author: ychuang
"""

import os
import sys

import flux_plot as fp
import B2pl_method as bm





def main(plotall):

    print("Initializing SOLPSxport")
    xp = fp.fluxplot(workdir = os.getcwd())
    
    print("Running getSOLPSlast10Profs")
    xp.getSOLPSlast10Profs(plotit = False, use_existing_last10 = False)
    
    print("Getting flux profiles")
    xp.getSOLPSfluxProfs(plotit = plotall)


# --- Launch main() ----------------------------------------------------------------------


if __name__ == '__main__':
    import argparse

    py3_9 = (sys.version_info[0] >= 3 and sys.version_info[1] >= 9)

    parser = argparse.ArgumentParser(description='Generate new b2.transport.inputfile files for SOLPS',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--plotall', help='plot flux profile', type= bool, default=None)
    
    args = parser.parse_args()
    
    main(plotall= args.plotall)
    