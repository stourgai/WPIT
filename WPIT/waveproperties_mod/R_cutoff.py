import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

#####waveproperties_mod.R_cutoff##############################

#Description:Cut-off frequency for R-mode waves
#Inputs:
# wpe: electron plasma frequency
# wce: electron cyclotron frequency
#Outputs:
# tmp: Cut-off frequency for R-mode waves

##############################################################

def R_cutoff(wce,wpe):
    sqrfac=np.sqrt(wce*wce+4*wpe*wpe)
    tmp=0.5*(wce+sqrfac)
    return tmp