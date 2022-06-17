import numpy as np


#####parallel_EMIC.nonlinear_H###############################################

#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# ppar_arg:parallel to the geomagnetic field momentum
# kappa_arg: parallel to the geomagnetic field wave number
# gamma_arg: Lorentz factor
# m_arg: particle mass
# wce_arg: electron gyrofrequency
# dkp_dt_arg: time derivative of the parallel compoment of the wave number
# dwcdz_arg: spatial derivative of the gyrofrequency along the field line
# dwdt_arg:time derivative of the wave frequency

#Outputs:
# tmp: H

#############################################################################

def nonlinear_H(pper_arg,ppar_arg,kappa_arg,gamma_arg,m_arg,wce_arg,dk_dt_arg,dwcdz_arg,dwdt_arg):
    dwc_dt=(ppar_arg/(gamma_arg*m_arg))*dwcdz_arg
    fac1=-(1/gamma_arg)*dwc_dt
    fac2=(ppar_arg/(gamma_arg*m_arg))*dk_dt_arg
    fac3=-((kappa_arg*pper_arg*pper_arg)/(2*gamma_arg*gamma_arg*m_arg*m_arg*wce_arg))*dwcdz_arg
    tmp=fac1+fac2+fac3-dwdt_arg
    return tmp