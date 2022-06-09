import numpy as np
from scipy.special import jn 

from environment_mod import const

def daeqdt(ppar_arg,pper_arg,eta_arg,gamma_arg,w_wave_arg,kappa_arg,Bw_arg,aeq_arg,alpha_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    aeq_rk=((q_arg*Bw_arg)/(p_mag**2))*(np.tan(aeq_arg)/np.tan(alpha_arg))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-(pper_arg**2/(gamma_arg*m_arg)))*np.sin(eta_arg)
    return aeq_rk