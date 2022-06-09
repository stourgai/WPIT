import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

def wave_amplitudes_bell(mu_arg,P_arg,D_arg,S_arg,Byw_arg,theta_arg):

    mu_sq_arg=mu_arg**2
    fac1= (P_arg-mu_sq_arg*(np.sin(theta_arg)**2)) 
    Byw_arg=Byw_arg

    Bxw_arg=(-(D_arg*fac1)/(P_arg*(S_arg-mu_arg**2)))*Byw_arg
    Bzw_arg=((D_arg*np.sin(theta_arg)*fac1)/(P_arg*np.cos(theta_arg)*(S_arg-mu_arg**2)))*Byw_arg
    Exw_arg=((const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg))*Byw_arg)
    Eyw_arg=((D_arg*const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg)*(mu_arg**2-S_arg)))*Byw_arg
    Ezw_arg=(-(const.c_light*mu_arg*np.sin(theta_arg))/P_arg)*Byw_arg
    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg