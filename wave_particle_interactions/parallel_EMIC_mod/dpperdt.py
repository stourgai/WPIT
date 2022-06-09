import numpy as np
from scipy.special import jn 

from environment_mod import const

def dpperdt(ppar_arg,pper_arg,eta_arg,Bw_arg,gamma_arg,wmega_arg,kappa_arg,wc_arg,dwcds_arg,q_arg,m_arg):
    mrk1=(q_arg*Bw_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.sin(eta_arg)
    mrk2=((pper_arg*ppar_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    mrk=mrk1+mrk2
    return mrk