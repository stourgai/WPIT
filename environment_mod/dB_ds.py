import numpy as np
from environment_mod import const

#####environment_mod.dB_ds############################################

#Description:Routine to calculate the gradient of magnetic field strength with respect to the distance along the magnetic field line
#Inputs:
# B_arg: magnetic field in T
# lamda_arg: magnetic latitude in rad
# L_arg: L shell
#Outputs:
# dB_ds_arg: gradient of magnetic field along the field line

#############################################################################

def dB_ds(B_arg,lamda_arg,L_arg):

    slat = np.sin(lamda_arg)
    clat = np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    dB_ds_arg = (3.0*B_arg/(L_arg*const.Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))

    return dB_ds_arg