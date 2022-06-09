import numpy as np
from scipy.special import jn 

from environment_mod import const

def dppardt(pper_arg,eta_arg,wc_arg,Bw_arg,gamma_arg,dwcds_arg,q_arg,m_arg):
    lrk1=((q_arg*Bw_arg)/(gamma_arg*m_arg))*pper_arg*np.sin(eta_arg)
    lrk2=((pper_arg*pper_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    lrk=lrk1-lrk2
    return lrk
