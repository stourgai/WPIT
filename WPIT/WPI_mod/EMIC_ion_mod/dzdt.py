import numpy as np

#####EMIC_ion_mod.dzdt###############################################

#Description:Routine to calculate the time derivative of the distance along the field line z
#Inputs:
# gamma_arg: Lorentz factor
# ppar_arg: parallel to the geomagnetic field momentum
# mi_arg: particle mass

#Outputs:
# tmp: 𝑑z𝑑𝑡

#############################################################################

def dzdt(ppar_arg,gamma_arg,mi_arg):
    tmp=ppar_arg/(mi_arg*gamma_arg)
    return tmp