import numpy as np


#####parallel_EMIC.dalphadt###############################################

#Description:Routine to calculate the time derivative of the local pitch angle
#Inputs:
# ppar_arg: parallel to the geomagnetic field momentum
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave particle phase
# gamma_arg: Lorentz factor
# wc_arg: particle cyclotron frequency
# dwcds_arg: derivative of cyclotron frequency with respect to the distance along the magnetic field line
# wmega_wave_arg: wave frequency in rad/s
# kappa_arg: wave number
# Bw_arg: wave magnetic field in T
# m_arg: particle mass
# q_arg: particle charge

#Outputs:
# prk: 𝑑𝛼𝑑𝑡

#############################################################################

def dalphadt(ppar_arg,pper_arg,eta_arg,gamma_arg,wce_arg,dwcds_arg,w_wave_arg,kappa_arg,Bw_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    a_rk1=((q_arg*Bw_arg)/(p_mag**2))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-((pper_arg**2)/(gamma_arg*m_arg)))*np.sin(eta_arg)
    a_rk2=(pper_arg/(2*gamma_arg*m_arg*wce_arg))*dwcds_arg
    a_rk=a_rk1+a_rk2
    return a_rk