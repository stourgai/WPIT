import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np

def L_cutoff(wce,wpe):
    sqrfac=np.sqrt(wce*wce+4*wpe*wpe)
    tmp=0.5*(-wce+sqrfac)
    return tmp