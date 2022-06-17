import numpy as np


#####parallel_EMIC.daeqdt###############################################

#Description:Routine to calculate the time derivative of the equatorial pitch angle
#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# ppar_arg: parallel to the geomagnetic field momentum
# eta_arg: wave-particle phase in rad
# gamma_arg: Lorentz factor
# wmega_wave_arg: wave frequency in rad/s
# kappa_arg: wave number
# Bw_arg: wave magnetic field in T
# aeq_arg: equatorial pitch angle
# alpha_arg: local pitch angle
# m_arg: particle mass
# q_arg: particle charge

#Outputs:
# tmp: 𝑑𝛼𝑒𝑞𝑑𝑡

#############################################################################


def daeqdt(ppar_arg,pper_arg,eta_arg,gamma_arg,w_wave_arg,kappa_arg,Bw_arg,aeq_arg,alpha_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    aeq_rk=((q_arg*Bw_arg)/(p_mag**2))*(np.tan(aeq_arg)/np.tan(alpha_arg))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-(pper_arg**2/(gamma_arg*m_arg)))*np.sin(eta_arg)
    return aeq_rk