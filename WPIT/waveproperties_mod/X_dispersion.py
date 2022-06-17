import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np

#####waveproperties_mod.X_dispersion##############################

#Description:Dispersion relation of X-mode (extra-ordinary) wave
#Inputs:
# w: wave frequency
# wpe: electron plasma frequency
# wlh: lower hybrid resonance frequency
#Outputs:
# nsq_tmp: squared refractive index
# kappa_tmp: wave number
#####################################################################

def X_dispersion(w,wpe,wlh):
    fac1=((wpe*wpe)/(w*w))
    fac2=((w*w-wpe*wpe)/(w*w-wlh*wlh))
    nsq_tmp=1-fac1*fac2




    tmp1=w*w/(const.c_light*const.c_light)
    tmpsq=tmp1*nsq_tmp
    kappa_tmp=np.sqrt(tmpsq)
    return nsq_tmp,kappa_tmp