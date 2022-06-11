import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

def warm_stix_parameters(S_cold,D_cold,P_cold,Te,Ti,mu_warm,Ke,KH,KHe,KO):
    Te_kelvin=Te*11604.5250061598  
    Ti_kelvin=Ti*11604.5250061598  
    qeT=(const.kb*Te_kelvin)/(const.me*const.c_light*const.c_light)
    qHT=(const.kb*Ti_kelvin)/(const.mH*const.c_light*const.c_light)
    qHeT=(const.kb*Ti_kelvin)/(const.mHe*const.c_light*const.c_light)
    qOT=(const.kb*Ti_kelvin)/(const.mO*const.c_light*const.c_light)    
    
    taue=qeT*mu_warm*mu_warm
    tauH=qHT*mu_warm*mu_warm
    tauHe=qHeT*mu_warm*mu_warm
    tauO=qOT*mu_warm*mu_warm
    
    K11e=Ke[0]
    K12e=Ke[1]
    K33e=Ke[8]
    K11H=KH[0]
    K11He=KHe[0]
    K11O=KO[0]
    
    S_warm=S_cold+taue*K11e+tauH*K11H+tauHe*K11He+tauO*K11O
    D_warm=D_cold+taue*K12e+tauH*K11H+tauHe*K11He+tauO*K11O
    P_warm=P_cold+taue*K33e+tauH*K11H+tauHe*K11He+tauO*K11O
    R_warm=S_warm+D_warm
    L_warm=S_warm-D_warm
    
    return S_warm,D_warm,P_warm,R_warm,L_warm