import numpy as np

#####parallel_EMIC.detadt###############################################

#Description:Routine to calculate the time derivative of the wave particle phase 𝜂
#Inputs:
# ppar_arg: parallel to the geomagnetic field momentum
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave particle phase
# Bw_arg: wave magnetic field in T
# wmega_arg: wave frequency in rad/s
# kappa_arg: wave number
# wc_arg: particle cyclotron frequency
# gamma_arg: Lorentz factor
# q_arg: particle charge
# m_arg: particle mass

#Outputs:
# mrk: 𝑑𝜂𝑑𝑡

#############################################################################

def detadt(ppar_arg,pper_arg,eta_arg,Bw_arg,wmega_arg,kappa_arg,wc_arg,gamma_arg,q_arg,m_arg,s):
    nrk1=((q_arg*Bw_arg)/pper_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.cos(eta_arg)
    nrk2a=(kappa_arg*ppar_arg)/(gamma_arg*m_arg)
    nrk2b=-wmega_arg
    if s=='ion':
        fac=1
    if s=='electron':
        fac=-1
    nrk2c=fac*(wc_arg/gamma_arg)
    nrk=nrk1+nrk2a+nrk2b+nrk2c
    return nrk