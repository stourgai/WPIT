import numpy as np
from environment_mod import const
def omega_plasma(n_arg,q_arg,m_arg):
    omegap_tmp=np.sqrt((n_arg*q_arg*q_arg)/(m_arg*const.epsilon0))

    return omegap_tmp