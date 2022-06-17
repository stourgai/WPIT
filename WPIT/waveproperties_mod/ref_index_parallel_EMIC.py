import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


import numpy as np

def ref_index_parallel_EMIC(wmega_wave,wpe,wce,wpH,wcH,wpO,wcO,wpHe,wcHe):
    if wpO==0 or wcO==0:
        Rfac_O=0
    else:
        Rfac_O=(wpO*wpO)/(wmega_wave*(wmega_wave-wcO))
    if wpHe==0 or wcHe==0:
        Rfac_He=0
    else:
        Rfac_He=(wpHe*wpHe)/(wmega_wave*(wmega_wave-wcHe)) 
        
    Rfac_e=(wpe*wpe)/(wmega_wave*(wmega_wave+wce))
    Rfac_H=(wpH*wpH)/(wmega_wave*(wmega_wave-wcH))
    

    R_tmp=1-Rfac_e-Rfac_H-Rfac_O-Rfac_He
    return R_tmp
