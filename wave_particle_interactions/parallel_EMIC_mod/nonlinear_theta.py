import numpy as np
from scipy.special import jn 

from environment_mod import const

def nonlinear_theta(pper_arg,ppar_arg,Bw_arg,kappa_arg,gamma_arg,m_arg,q_arg,wce_arg,w_arg):
    fac1=q_arg*Bw_arg*pper_arg
    fac2a=kappa_arg/(gamma_arg*gamma_arg*m_arg*m_arg)
    fac2b=(wce_arg-((kappa_arg*ppar_arg)/m_arg))*(w_arg/(kappa_arg*gamma_arg*gamma_arg*gamma_arg*m_arg*const.c_light*const.c_light))
    tmp=fac1*(fac2a+fac2b)
    tmpwtsq_arg=np.abs(tmp)
    return tmp,tmpwtsq_arg