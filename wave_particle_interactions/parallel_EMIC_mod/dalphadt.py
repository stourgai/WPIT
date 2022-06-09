import numpy as np
from scipy.special import jn 

from environment_mod import const

def dalphadt(ppar_arg,pper_arg,eta_arg,gamma_arg,wce_arg,dwcds_arg,w_wave_arg,kappa_arg,Bw_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    a_rk1=((q_arg*Bw_arg)/(p_mag**2))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-((pper_arg**2)/(gamma_arg*m_arg)))*np.sin(eta_arg)
    a_rk2=(pper_arg/(2*gamma_arg*m_arg*wce_arg))*dwcds_arg
    a_rk=a_rk1+a_rk2
    return a_rk