import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np

def light_dispersion(w,wpe,wce):
    nsq_tmp=(w*w-wpe*wpe)/(w*w)

    tmpsq=(w*w-wpe*wpe)/(const.c_light*const.c_light)
    kappa_tmp=np.sqrt(tmpsq)
    return nsq_tmp,kappa_tmp