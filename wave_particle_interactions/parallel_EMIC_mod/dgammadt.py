import numpy as np
from scipy.special import jn 

from environment_mod import const

def dgammadt(pper_arg,eta_arg,gamma_arg,Bw_arg,kappa_arg,w_arg,q_arg,m_arg):
    fac1=1/(gamma_arg*m_arg*m_arg*const.c_light*const.c_light)
    fac2=q_arg*pper_arg*Bw_arg*(w_arg/kappa_arg)*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp