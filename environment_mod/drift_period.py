import numpy as np
from environment_mod import const

#####environment_mod.drift_period############################################

#Description:Calculate the drift period of a trapped particle
#Inputs:
# pper_arg: paerpendicular to the magnetic field momentum component
# B_arg: Magnetic field strength
# m_arg: Particle mass
#Outputs:
# tmp: Drift period in s

#############################################################################

def drift_period(B_arg,m_arg,v_arg,L_arg,aeq_arg):
    fac1=(2*np.pi*const.qe*B_arg*((const.Re)**3))/(m_arg*v_arg*v_arg)
    fac2=1/(L_arg*const.Re)
    fac3=1-(1/3)*((np.sin(aeq_arg))**0.62)
    tmp=fac1*fac2*fac3
    return tmp
