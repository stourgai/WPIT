import numpy as np


#####parallel_EMIC.dpperdt###############################################

#Description:Routine to calculate the time derivative of the perpendicular momentum
#Inputs:
# ppar_arg: parallel to the geomagnetic field momentum
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave particle phase
# Bw_arg: wave magnetic field in T
# gamma_arg: Lorentz factor
# wmega_arg: wave frequency in rad/s
# kappa_arg: wave number
# wc_arg: particle cyclotron frequency
# dwcds_arg: derivative of cyclotron frequency with respect to the distance along the magnetic field line
# q_arg: particle charge
# m_arg: particle mass

#Outputs:
# mrk: 𝑑𝑝⊥𝑑𝑡

#############################################################################

def dpperdt(ppar_arg,pper_arg,eta_arg,Bw_arg,gamma_arg,wmega_arg,kappa_arg,wc_arg,dwcds_arg,q_arg,m_arg):
    mrk1=(q_arg*Bw_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.sin(eta_arg)
    mrk2=((pper_arg*ppar_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    mrk=mrk1+mrk2
    return mrk