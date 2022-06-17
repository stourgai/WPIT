import numpy as np


#####parallel_EMIC.dEkdt###############################################

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
# tmp: 𝑑𝐸𝑘𝑑𝑡

#############################################################################

def dEkdt(pper_arg,eta_arg,gamma_arg,Bw_arg,kappa_arg,w_arg,q_arg,m_arg):
    fac1=1/(gamma_arg*m_arg)
    fac2=q_arg*pper_arg*Bw_arg*(w_arg/kappa_arg)*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp