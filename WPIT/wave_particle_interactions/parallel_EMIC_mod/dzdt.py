import numpy as np
from scipy.special import jn 

from environment_mod import const

#####parallel_EMIC.dzdt###############################################

#Description:Routine to calculate the time derivative of the distance along the field line z
#Inputs:
# gamma_arg: Lorentz factor
# ppar_arg: parallel to the geomagnetic field momentum
# mi_arg: particle mass

#Outputs:
# tmp: 𝑑𝑧𝑑𝑡

#############################################################################


def dzdt(ppar_arg,gamma_arg,mi_arg):
    krk=ppar_arg/(gamma_arg*mi_arg)
    return krk
