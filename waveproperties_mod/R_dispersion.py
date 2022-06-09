import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np

def R_dispersion(w,wpe,wce,wpH,wcH,wpHe,wcHe,wpO,wcO):
    nsq_tmpe=(wpe*wpe)/(w*(w-wce))
    nsq_tmpH=(wpH*wpH)/(w*(w+wcH))
    nsq_tmpHe=(wpHe*wpHe)/(w*(w+wcHe))
    nsq_tmpO=(wpO*wpO)/(w*(w+wcO))

    nsq_tmp=1-nsq_tmpe-nsq_tmpH-nsq_tmpHe-nsq_tmpO
        
    tmp1=w*w/(const.c_light*const.c_light)
    tmpsq=tmp1*nsq_tmp
    kappa_tmp=np.sqrt(tmpsq)
    return nsq_tmp,kappa_tmp