import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

#####waveproperties_mod.wave_amplitudes_li##############

#Description:Routine to calculate the wave amplitudes
#Inputs:
# mu: the refractive index
# P: Stix P parameter
# D: Stix D parameter
# S: Stix S parameter
# Bw_tot_li: magnitude of the magentic field of the wave in T
# psi: wave normal angle in rad
#Outputs:
# Bxw_arg:Bx component of the wave in T
# Byw_arg:By component of the wave in T
# Bzw_arg:Bz component of the wave in T
# Exw_arg:Ex component of the wave in V/m
# Eyw_arg:Ey component of the wave in V/m
# Ezw_arg:Ez component of the wave in V/m
#################################################################

def wave_amplitudes_li(mu,P,D,S,Bw_tot_li,psi):
    #Li uses a different approach in defining the wave fields (see notebook for more details)

    I_w=(Bw_tot_li/(mu*np.sqrt(D*D*(P-mu*mu*np.sin(psi)*np.sin(psi))*(P-mu*mu*np.sin(psi)*np.sin(psi))
                           +P*P*np.cos(psi)*np.cos(psi)*(S-mu*mu)*(S-mu*mu))))
    mu_sq_li=mu*mu
    fac1= (P-mu_sq_li*np.sin(psi)*np.sin(psi)) 
    fac2=S-mu*mu
    Exw_li=const.c_light*I_w*fac1*fac2
    Eyw_li=const.c_light*I_w*D*fac1
    Ezw_li=-const.c_light*I_w*mu_sq_li*np.cos(psi)*np.sin(psi)*fac2
    Bxw_li=-I_w*D*np.cos(psi)*fac1*mu
    Byw_li=I_w*P*np.cos(psi)*fac2*mu
    Bzw_li=I_w*D*np.sin(psi)*fac1*mu
                
    return Bxw_li, Byw_li, Bzw_li, Exw_li, Eyw_li, Ezw_li
