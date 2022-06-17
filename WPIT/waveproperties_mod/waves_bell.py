import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

def waves_bell(Bxwc_arg, Bywc_arg, Bzwc_arg, Exwc_arg, Eywc_arg, Ezwc_arg,time_arg,kz_arg,zeta_arg,kx_arg,chi_arg,w_wave_arg):
    Phi_arg=(w_wave_arg*time_arg-kz_arg*zeta_arg-kx_arg*chi_arg)   #should i add kx*x ?

    Exw_arg=-Exwc_arg*np.sin(Phi_arg)  
    Eyw_arg=Eywc_arg*np.cos(Phi_arg)
    Ezw_arg=-Ezwc_arg*np.sin(Phi_arg)
    Bxw_arg=Bxwc_arg*np.cos(Phi_arg)  
    Byw_arg=Bywc_arg*np.sin(Phi_arg)
    Bzw_arg=-Bzwc_arg*np.cos(Phi_arg)

    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg, Phi_arg
