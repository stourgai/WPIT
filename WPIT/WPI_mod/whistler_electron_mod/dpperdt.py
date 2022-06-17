import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def dpperdt(ppar_arg,pper_arg,eta_arg,w1_arg,w2_arg,beta_arg,gamma_arg,R1_arg,R2_arg,m_res_arg,wce_arg,dwds_arg):
    mrk=(-pow((-1),(m_res_arg-1))*(w1_arg*((ppar_arg/gamma_arg)+(const.me*R1_arg))*jn((m_res_arg-1),beta_arg)-
            w2_arg*((ppar_arg/gamma_arg)-(const.me*R2_arg))*jn((m_res_arg+1),beta_arg))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg*ppar_arg)/(2*wce_arg))*dwds_arg)
    return mrk