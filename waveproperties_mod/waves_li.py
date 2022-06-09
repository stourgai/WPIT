import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

def waves_li(Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc,time,kz,kx,zeta,chi,w_wave):
    Phi=(w_wave*time-kz*zeta-kx*chi) 

    Exw=Exwc*np.sin(Phi)  
    Eyw=Eywc*np.cos(Phi)
    Ezw=Ezwc*np.sin(Phi)
    Bxw=Bxwc*np.cos(Phi)  
    Byw=Bywc*np.sin(Phi)
    Bzw=Bzwc*np.cos(Phi)

    return Bxw, Byw, Bzw, Exw, Eyw, Ezw, Phi
