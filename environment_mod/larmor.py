import numpy as np
from environment_mod import const

#####environment_mod.larmor############################################

#Description:Routine to calculatethe larmor radius
#Inputs:
# uperp_arg: perpendicular velocity in m/s
# gamma_arg: Lorentz factor
# B_arg: Magnetic field in T
# ms_arg: particle mass
# qs_arg: particle charge
#Outputs:
# rl_tmp: Larmor radius in m

#############################################################################

def larmor(uperp_arg,gamma_arg,B_arg,ms_arg,qs_arg):
    rl_tmp=gamma_arg*ms_arg*uperp_arg/(qs_arg*B_arg)
    return rl_tmp