import numpy as np
from environment_mod import const

#####environment_mod.omega_plasma############################################

#Description:Routine to calculate the plasma frequency
#Inputs:
# n_arg: particle number density in 𝑚−3
# q_arg: particle charge in Cb
# m_arg: particle mass in kg
#Outputs:
# omegap_tmp: plasma frequency in rad/s

#############################################################################

def omega_plasma(n_arg,q_arg,m_arg):
    omegap_tmp=np.sqrt((n_arg*q_arg*q_arg)/(m_arg*const.epsilon0))

    return omegap_tmp