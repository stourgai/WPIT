"""
waveproperties_mod.ref_index_parallel_EMIC

**Description**:
_____________________________________________________________________________________________________________________

Calculate the refractive index of parallel propagating EMIC waves
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

wmega_wave: wave frequency in rad/s

wpe: electron plasma frequency in rad/s

wce: electron cyclotron frequency in rad/s

wpH: H+ plasma frequency in rad /s

wcH: H+ cyclotron frequency in rad/s

wpO: O+ plasma frequency in rad /s

wcO: O+ cyclotron frequency in rad/s

wpHe: He+ plasma frequency in rad /s

wcHe: He+ cyclotron frequency in rad/s

______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

R_tmp: refractive index
______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Summers, D. and Thorne, R. M. (2003). Relativistic electron pitch-angle scattering by electromagnetic ion
cyclotron waves during geomagnetic storms. Journal of Geophysical Research: Space Physics 108
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np

def refr_index_parallel_EMIC(wmega_wave,wpe,wce,wpH,wcH,wpO,wcO,wpHe,wcHe):
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
