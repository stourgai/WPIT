import numpy as np
from scipy.special import jn 


def nonlinear_theta(C0,Cp1,Cm1,m_res,beta):
    tmp=C0*jn((m_res),beta)+Cp1*jn((m_res+1),beta)+Cm1*jn((m_res-1),beta)
    tmpwtsq=np.abs(tmp)
    return tmp,tmpwtsq