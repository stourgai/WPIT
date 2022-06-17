import numpy as np

from WPIT.Environment_mod import const

#####parallel_EMIC.nonlinear_theta###############################################

#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# ppar_arg: parallel to the geomagnetic field momentum
# Bw_arg: wave magnetic field in T
# kappa_arg: wave number
# gamma_arg: Lorentz factor
# m_arg: particle mass
# q_arg: particle charge
# wce_arg: cycotron frequency in rad/s
# w_arg: wave frequency in rad/s

#Outputs:
# tmp: non linear parameter 𝜃
# tmpwtsq_arg: 𝜔2𝑡

#############################################################################


def nonlinear_theta(pper_arg,ppar_arg,Bw_arg,kappa_arg,gamma_arg,m_arg,q_arg,wce_arg,w_arg):
    fac1=q_arg*Bw_arg*pper_arg
    fac2a=kappa_arg/(gamma_arg*gamma_arg*m_arg*m_arg)
    fac2b=(wce_arg-((kappa_arg*ppar_arg)/m_arg))*(w_arg/(kappa_arg*gamma_arg*gamma_arg*gamma_arg*m_arg*const.c_light*const.c_light))
    tmp=fac1*(fac2a+fac2b)
    tmpwtsq_arg=np.abs(tmp)
    return tmp,tmpwtsq_arg