import numpy as np
from scipy.special import jn 

from environment_mod import const

def detadt(ppar_arg,pper_arg,eta_arg,Bw_arg,wmega_arg,kappa_arg,wc_arg,gamma_arg,q_arg,m_arg):
    nrk1=((q_arg*Bw_arg)/pper_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.cos(eta_arg)
    nrk2a=(kappa_arg*ppar_arg)/(gamma_arg*m_arg)
    nrk2b=-wmega_arg
    nrk2c=-wc_arg/gamma_arg
    nrk=nrk1+nrk2a+nrk2b+nrk2c
    return nrk