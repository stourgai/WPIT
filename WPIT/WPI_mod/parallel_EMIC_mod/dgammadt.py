import numpy as np


from WPIT.Environment_mod import const

#####parallel_EMIC.dgammadt###############################################

#Description:Routine to calculate the time derivative of the equatorial pitch angle
#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave-particle phase in rad
# gamma_arg: Lorentz factor
# Bw_arg: wave magnetic field in T
# kappa_arg: wave number
# wmega_wave_arg: wave frequency in rad/s
# q_arg: particle charge
# m_arg: particle mass

#Outputs:
# tmp: 𝑑𝛾𝑑𝑡

#############################################################################

def dgammadt(pper_arg,eta_arg,gamma_arg,Bw_arg,kappa_arg,w_arg,q_arg,m_arg):
    fac1=1/(gamma_arg*m_arg*m_arg*const.c_light*const.c_light)
    fac2=q_arg*pper_arg*Bw_arg*(w_arg/kappa_arg)*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp