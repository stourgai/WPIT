import numpy as np

#####environment_mod.Bmag_dipole###############################################

#Description:Routine to calculate the geomagnetic dipole field strength
#Inputs:
# L_arg: L shell
# lamda_arg: magnetic latitude in rad
#Outputs:
# Bmag: Geomagnetic field in T

#############################################################################

def Bmag_dipole(L_arg,lamda_arg):
    
   
    B0=31200*10**(-9)
    slat=np.sin(lamda_arg)
    clat=np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    Bmag=(B0/(L_arg**3))*slat_term/(clat**6)

    return Bmag
