import numpy as np

#####environment_mod.mu_adiabatic############################################

#Description:Routine to calculatethe the first adiabatic invariant
#Inputs:
# pper_arg: perpendicualr to the magnetic field component of momentum
# B_arg: magnetic field strength
# m_arg: particle mass
#Outputs:
# tmp: first adiabatic invariant

#############################################################################

def mu_adiabatic(pper_arg,B_arg,m_arg):
    tmp=(pper_arg*pper_arg)/(2*m_arg*B_arg)
    return tmp
