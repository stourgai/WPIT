"""
waveproperties_mod.warm_dielectric_tensor_warm

**Description**:
_____________________________________________________________________________________________________________________

Calculate the dielectric tensor with warm plasma corrections
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

psi_arg: wave normal angle in rad
_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

K_e:the electron warm dielectric tensor compoments

K_H:the hydrogen warm dielectric tensor compoments

K_He:the helium warm dielectric tensor compoments

K_O: the oxygen warm dielectric tensor compoments 
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Maxworth, A. and Gołkowski, M. (2017). Magnetospheric whistler mode ray tracing in a warm background
plasma with finite electron and ion temperature. Journal of Geophysical Research: Space Physics 122,
7323–7335
_________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np


from WPIT.Environment_mod import const


def warm_dielectric_tensor(w_wave,wce_arg,wpe_arg,wcH_arg,wpH_arg,wcHe_arg,wpHe_arg,wcO_arg,wpO_arg,psi_arg):
    Xe=(wpe_arg*wpe_arg)/(w_wave*w_wave)
    XH=(wpH_arg*wpH_arg)/(w_wave*w_wave)    
    XHe=(wpHe_arg*wpHe_arg)/(w_wave*w_wave)    
    XO=(wpO_arg*wpO_arg)/(w_wave*w_wave)  
    
    Ye=-wce_arg/w_wave
    YH=wcH_arg/w_wave
    YHe=wcHe_arg/w_wave
    YO=wcO_arg/w_wave
    
    K_11_e_fac1=-Xe/(1-Ye*Ye)
    K_11_e_fac2=(3*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*Ye*Ye)
    K_11_e_fac3=((1+3*Ye*Ye)*np.cos(psi_arg)*np.cos(psi_arg))/((1-Ye*Ye)*(1-Ye*Ye))
    
    K_11_e=K_11_e_fac1*(K_11_e_fac2+K_11_e_fac3)
    
    K_12_e_fac1=1j*Xe/(1-Ye*Ye)
    K_12_e_fac2=(6*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*Ye*Ye)
    K_12_e_fac3=((3+Ye*Ye)*np.cos(psi_arg)*np.cos(psi_arg))/((1-Ye*Ye)*(1-Ye*Ye))
  
    K_12_e=K_12_e_fac1*(K_12_e_fac2+K_12_e_fac3)
    
    K_13_e=(2*Xe*Ye*np.sin(psi_arg)*np.cos(psi_arg))/((1-Ye*Ye)*(1-Ye*Ye))
    
    K_21_e=-K_12_e
    
    K_22_e_fac1=-Xe/(1-Ye*Ye)
    K_22_e_fac2=((1+8*Ye*Ye)*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*Ye*Ye)
    K_22_e_fac3=((1+3*Ye*Ye)*np.cos(psi_arg)*np.cos(psi_arg))/((1-Ye*Ye)*(1-Ye*Ye))
    
    K_22_e=K_22_e_fac1*(K_22_e_fac2+K_22_e_fac3)
    
    K_23_e=(-1j*Xe*Ye*(3-Ye*Ye)*np.cos(psi_arg)*np.sin(psi_arg))/((1-Ye*Ye)*(1-Ye*Ye))
    
    K_31_e=K_13_e
    
    K_32_e=-K_23_e
    
    K_33_e_fac1=-Xe
    K_33_e_fac2=3*np.cos(psi_arg)*np.cos(psi_arg)
    K_33_e_fac3=(np.sin(psi_arg)*np.sin(psi_arg)/(1-Ye*Ye))
    
    K_33_e=K_33_e_fac1*(K_33_e_fac2+K_33_e_fac3)
    
    K_11_H_fac1=-XH/(1-YH*YH)
    K_11_H_fac2=(3*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YH*YH)
    K_11_H_fac3=((1+3*YH*YH)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YH*YH)*(1-YH*YH))
    
    K_11_H=K_11_H_fac1*(K_11_H_fac2+K_11_H_fac3)
    
    K_12_H_fac1=1j*XH/(1-YH*YH)
    K_12_H_fac2=(6*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YH*YH)
    K_12_H_fac3=((3+YH*YH)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YH*YH)*(1-YH*YH))
  
    K_12_H=K_12_H_fac1*(K_12_H_fac2+K_12_H_fac3)
    
    K_13_H=(2*XH*YH*np.sin(psi_arg)*np.cos(psi_arg))/((1-YH*YH)*(1-YH*YH))
    
    K_21_H=-K_12_H
    
    K_22_H_fac1=-XH/(1-YH*YH)
    K_22_H_fac2=((1+8*YH*YH)*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YH*YH)
    K_22_H_fac3=((1+3*YH*YH)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YH*YH)*(1-YH*YH))
    
    K_22_H=K_22_H_fac1*(K_22_H_fac2+K_22_H_fac3)
    
    K_23_H=(-1j*XH*YH*(3-YH*YH)*np.cos(psi_arg)*np.sin(psi_arg))/((1-YH*YH)*(1-YH*YH))
    
    K_31_H=K_13_H
    
    K_32_H=-K_23_H
    
    K_33_H_fac1=-XH
    K_33_H_fac2=3*np.cos(psi_arg)*np.cos(psi_arg)
    K_33_H_fac3=(np.sin(psi_arg)*np.sin(psi_arg)/(1-YH*YH))
    
    K_33_H=K_33_H_fac1*(K_33_H_fac2+K_33_H_fac3)
    
    K_11_He_fac1=-XHe/(1-YHe*YHe)
    K_11_He_fac2=(3*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YHe*YHe)
    K_11_He_fac3=((1+3*YHe*YHe)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YHe*YHe)*(1-YHe*YHe))
    
    K_11_He=K_11_He_fac1*(K_11_He_fac2+K_11_He_fac3)
    
    K_12_He_fac1=1j*XHe/(1-YHe*YHe)
    K_12_He_fac2=(6*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YHe*YHe)
    K_12_He_fac3=((3+YHe*YHe)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YHe*YHe)*(1-YHe*YHe))
  
    K_12_He=K_12_He_fac1*(K_12_He_fac2+K_12_He_fac3)
    
    K_13_He=(2*XHe*YHe*np.sin(psi_arg)*np.cos(psi_arg))/((1-YHe*YHe)*(1-YHe*YHe))
    
    K_21_He=-K_12_He
    
    K_22_He_fac1=-XHe/(1-YHe*YHe)
    K_22_He_fac2=((1+8*YHe*YHe)*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YHe*YHe)
    K_22_He_fac3=((1+3*YHe*YHe)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YHe*YHe)*(1-YHe*YHe))
    
    K_22_He=K_22_He_fac1*(K_22_He_fac2+K_22_He_fac3)
    
    K_23_He=(-1j*XHe*YHe*(3-YHe*YHe)*np.cos(psi_arg)*np.sin(psi_arg))/((1-YHe*YHe)*(1-YHe*YHe))
    
    K_31_He=K_13_He
    
    K_32_He=-K_23_He
    
    K_33_He_fac1=-XHe
    K_33_He_fac2=3*np.cos(psi_arg)*np.cos(psi_arg)
    K_33_He_fac3=(np.sin(psi_arg)*np.sin(psi_arg)/(1-YHe*YHe))
    
    K_33_He=K_33_He_fac1*(K_33_He_fac2+K_33_He_fac3)
    
    K_11_O_fac1=-XO/(1-YO*YO)
    K_11_O_fac2=(3*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YO*YO)
    K_11_O_fac3=((1+3*YO*YO)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YO*YO)*(1-YO*YO))
    
    K_11_O=K_11_O_fac1*(K_11_O_fac2+K_11_O_fac3)
    
    K_12_O_fac1=1j*XO/(1-YO*YO)
    K_12_O_fac2=(6*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YO*YO)
    K_12_O_fac3=((3+YO*YO)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YO*YO)*(1-YO*YO))
  
    K_12_O=K_12_O_fac1*(K_12_O_fac2+K_12_O_fac3)
    
    K_13_O=(2*XO*YO*np.sin(psi_arg)*np.cos(psi_arg))/((1-YO*YO)*(1-YO*YO))
    
    K_21_O=-K_12_O
    
    K_22_O_fac1=-XO/(1-YO*YO)
    K_22_O_fac2=((1+8*YO*YO)*np.sin(psi_arg)*np.sin(psi_arg))/(1-4*YO*YO)
    K_22_O_fac3=((1+3*YO*YO)*np.cos(psi_arg)*np.cos(psi_arg))/((1-YO*YO)*(1-YO*YO))
    
    K_22_O=K_22_O_fac1*(K_22_O_fac2+K_22_O_fac3)
    
    K_23_O=(-1j*XO*YO*(3-YO*YO)*np.cos(psi_arg)*np.sin(psi_arg))/((1-YO*YO)*(1-YO*YO))
    
    K_31_O=K_13_O
    
    K_32_O=-K_23_O
    
    K_33_O_fac1=-XO
    K_33_O_fac2=3*np.cos(psi_arg)*np.cos(psi_arg)
    K_33_O_fac3=(np.sin(psi_arg)*np.sin(psi_arg)/(1-YO*YO))
    
    K_33_O=K_33_O_fac1*(K_33_O_fac2+K_33_O_fac3)
    
    K_11=K_11_e+K_11_H+K_11_He+K_11_O
    K_12=K_12_e+K_12_H+K_12_He+K_12_O
    K_13=K_13_e+K_13_H+K_13_He+K_13_O
    K_21=K_21_e+K_21_H+K_21_He+K_21_O
    K_22=K_22_e+K_22_H+K_22_He+K_22_O
    K_23=K_23_e+K_23_H+K_23_He+K_23_O
    K_31=K_31_e+K_31_H+K_31_He+K_31_O
    K_32=K_32_e+K_32_H+K_32_He+K_32_O
    K_33=K_33_e+K_33_H+K_33_He+K_33_O
    
    K_e=[K_11_e,K_12_e,K_13_e,K_21_e,K_22_e,K_23_e,K_31_e,K_32_e,K_33_e]
    K_H=[K_11_H,K_12_H,K_13_H,K_21_H,K_22_H,K_23_H,K_31_H,K_32_H,K_33_H]
    K_He=[K_11_He,K_12_He,K_13_He,K_21_He,K_22_He,K_23_He,K_31_He,K_32_He,K_33_He]
    K_O=[K_11_O,K_12_O,K_13_O,K_21_O,K_22_O,K_23_O,K_31_O,K_32_O,K_33_O]
    
    return K_e,K_H,K_He,K_O
