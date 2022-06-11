import numpy as np
from environment_mod import const

#####environment_mod.loss_cone_v2############################################

#Description:Routine to calculate the loss cone angle
#Inputs:
# **L_arg**: L shell
#Outputs:
# **tmp**: loss cone angle in rad

#############################################################################

def loss_cone_v2(L_arg,h_arg):
    zeta_m=(const.Re+h_arg)/(L_arg*const.Re)
    facsqrt=np.sqrt(1+3*(1-zeta_m))
    sinalc_sq=(zeta_m**3)/facsqrt
    sinalc=np.sqrt(sinalc_sq)
    tmp=np.arcsin(sinalc)

    return tmp
