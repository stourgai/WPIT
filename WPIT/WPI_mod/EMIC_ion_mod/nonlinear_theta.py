import numpy as np
from scipy.special import jn 


def nonlinear_theta(C0_arg,Cp1_arg,Cm1_arg,m_res_arg,beta_arg):
    tmp=C0_arg*jn((m_res_arg),beta_arg)+Cp1_arg*jn((m_res_arg+1),beta_arg)+Cm1_arg*jn((m_res_arg-1),beta_arg)
    tmpwtsq=np.abs(tmp)
    return tmp,tmpwtsq