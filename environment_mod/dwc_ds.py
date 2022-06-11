import numpy as np
from environment_mod import const

#####environment_mod.dwc_ds############################################

#Description:Routine to calculate the gradient of gyrofrequency with respect to the distance along the magnetic field line
#Inputs:
# wc_arg: gyrofrequency in rad/s
# lamda_arg: magnetic latitude in rad
# L_arg: L shell
#Outputs:
# dwce_ds_arg: gradient of gyrofrequency along the field line

#############################################################################

def dwc_ds(wc_arg,lamda_arg,L_arg):

    slat = np.sin(lamda_arg)
    clat = np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    dwce_ds_arg = (3.0*wc_arg/(L_arg*const.Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))
    return dwce_ds_arg