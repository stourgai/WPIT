import numpy as np


#####parallel_EMIC.dppardt###############################################

#Description:Routine to calculate the time derivative of the parallel momentum
#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave-particle phase in rad
# wc_arg: particle cyclotron frequency
# Bw_arg: wave magnetic field in T
# gamma_arg: Lorentz factor
# dwcds_arg: derivative of cyclotron frequency with respect to the distance along the magnetic field line
# q_arg: particle charge
# m_arg: particle mass

#Outputs:
# lrk: 𝑑𝑝∥𝑑𝑡

#############################################################################

def dppardt(pper_arg,eta_arg,wc_arg,Bw_arg,gamma_arg,dwcds_arg,q_arg,m_arg):
    lrk1=((q_arg*Bw_arg)/(gamma_arg*m_arg))*pper_arg*np.sin(eta_arg)
    lrk2=((pper_arg*pper_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    lrk=lrk1-lrk2
    return lrk
