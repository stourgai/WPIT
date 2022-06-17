import numpy as np

from WPIT.Environment_mod import const

#####parallel_EMIC.dlamdadt###############################################

#Description:Routine to calculate the time derivative of the magnetic latitude
#Inputs:
# ppar_arg: parallel to the geomagnetic field momentum
# gamma_arg: Lorentz factor
# L_arg: L shell

#Outputs:
# ork: 𝑑𝜆𝑑𝑡

#############################################################################

def dlamdadt(ppar_arg,lamda_arg,gamma_arg,m_arg,L_arg):
    ork=ppar_arg/(gamma_arg*m_arg*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork