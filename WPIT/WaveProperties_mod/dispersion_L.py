"""
waveproperties_mod.dispersion_L

**Description**:
_____________________________________________________________________________________________________________________

Dispersion relation of L-mode wave

_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w: wave frequency

wpe: electron plasma frequency

wce: electron cyclotron frequency

wpH: hydrogen plasma frequency

wcH: hydrogen cyclotron frequency

wpHe: helium plasma frequency (if available, else 0)

wcHe: helium cyclotron frequency (if available, else 0)

wpO: oxygen plasma frequency (if available, else 0)

wcO: oxygen cyclotron frequency (if available, else 0)

______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

nsq_tmp: squared refractive index

kappa_tmp: wave number
______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Swanson, D. G. (2012). Plasma waves (Elsevier)
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
from WPIT.Environment_mod import const


import numpy as np

def dispersion_L(w,wpe,wce,wpH,wcH,wpHe,wcHe,wpO,wcO):
    nsq_tmpe=(wpe*wpe)/(w*(w+wce))
    nsq_tmpH=(wpH*wpH)/(w*(w-wcH))
    nsq_tmpHe=(wpHe*wpHe)/(w*(w-wcHe))
    nsq_tmpO=(wpO*wpO)/(w*(w-wcO))

    nsq_tmp=1-nsq_tmpe-nsq_tmpH-nsq_tmpHe-nsq_tmpO    

    tmp1=w*w/(const.c_light*const.c_light)
    tmpsq=tmp1*nsq_tmp
    kappa_tmp=np.sqrt(tmpsq)

    return nsq_tmp,kappa_tmp
