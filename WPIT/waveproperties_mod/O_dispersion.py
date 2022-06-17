import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np


#####waveproperties_mod.O_dispersion##############################

#Description:Dispersion relation of O-mode (ordinary) wave
#Inputs:
# w: wave frequency
# wpe: electron plasma frequency
# wpH: hydrogen plasma frequency
# wpHe: helium plasma frequency (if available, else 0)
# wpO: oxygen plasma frequency (if available, else 0)
#Outputs:
# nsq_tmp: squared refractive index
# kappa_tmp: wave number
#####################################################################

def O_dispersion(w,wpe,wpH,wpHe,wpO):
    nsq_tmpe=(wpe*wpe)/(w*w)
    nsq_tmpH=(wpH*wpH)/(w*w)
    nsq_tmpHe=(wpHe*wpHe)/(w*w)
    nsq_tmpO=(wpO*wpO)/(w*w)

    nsq_tmp=1-nsq_tmpe-nsq_tmpH-nsq_tmpHe-nsq_tmpO    

    tmp1=w*w/(const.c_light*const.c_light)
    tmpsq=tmp1*nsq_tmp
    kappa_tmp=np.sqrt(tmpsq)
    
    return nsq_tmp,kappa_tmp