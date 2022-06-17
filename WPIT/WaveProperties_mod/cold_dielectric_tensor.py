"""
waveproperties_mod.dielectric_tensor_cold

**Description**:
_____________________________________________________________________________________________________________________

Calculate the dielectric tensor based on cold plasma theory
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w_wave: wave frequency

wce_arg: 𝑒 plasma frequency

wpe_arg: 𝑒 cyclotron frequency

wcH_arg: 𝐻+ cyclotron frequency

wpH_arg: 𝐻+ plasma frequency

wcHe_arg: 𝐻𝑒+ cyclotron frequency

wpHe_arg: 𝐻𝑒+ plasma frequency

wcO_arg: 𝑂+ cyclotron frequency

wpO_arg: 𝑂+ plasma frequency

_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

K_11,K_12,K_13,K_21,K_22,K_23,K_31,K_32,K_33: the cold dielectric tensor compoments 
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________
Maxworth, A. S., et al. "Raytracing study of source regions of whistler mode wave power distribution relative to the plasmapause." Journal of Geophysical Research: Space Physics 125.4 (2020): e2019JA027154.

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np

def cold_dielectric_tensor(w_wave,wce_arg,wpe_arg,wcH_arg,wpH_arg,wcHe_arg,wpHe_arg,wcO_arg,wpO_arg):
    X_e=(wpe_arg*wpe_arg)/(w_wave*w_wave)
    X_H=(wpH_arg*wpH_arg)/(w_wave*w_wave)    
    X_He=(wpHe_arg*wpHe_arg)/(w_wave*w_wave)    
    X_O=(wpO_arg*wpO_arg)/(w_wave*w_wave)  
    
    Y_e=-wce_arg/w_wave
    Y_H=wcH_arg/w_wave
    Y_He=wcHe_arg/w_wave
    Y_O=wcO_arg/w_wave
    
    K11_e=X_e/(1-Y_e*Y_e)
    K11_H=X_H/(1-Y_H*Y_H)
    K11_He=X_He/(1-Y_He*Y_He)
    K11_O=X_O/(1-Y_O*Y_O)
    
    K_11=1-K11_e-K11_H-K11_He-K11_O
    
    K12_e=1j*(X_e*Y_e)/(1-Y_e*Y_e)
    K12_H=1j*(X_H*Y_H)/(1-Y_H*Y_H)
    K12_He=1j*(X_He*Y_He)/(1-Y_He*Y_He)
    K12_O=1j*(X_O*Y_O)/(1-Y_O*Y_O)
    
    K_12=K12_e+K12_H+K12_He+K12_O    
    
    K_13=0
    
    K_21=-K_12
    
    K_22=K_11
    
    K_23=0
    
    K_31=0
    
    K_32=0
    
    K_33=1-X_e-X_H-X_He-X_O
    
    return K_11,K_12,K_13,K_21,K_22,K_23,K_31,K_32,K_33
