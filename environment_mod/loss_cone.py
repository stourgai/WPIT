import numpy as np

#####environment_mod.loss_cone###############################################

#Description:Routine to calculate the loss cone angle
#Inputs:
# L_arg: L shell
#Outputs:
# tmp: loss cone angle in rad

#############################################################################

def loss_cone(L_arg):
    fac=1/np.sqrt(L_arg**3*np.sqrt(4-3/L_arg))
    tmp=np.arcsin(fac)

    return tmp
