import numpy as np
from scipy.special import jn 

from environment_mod import const


def dppardt(pper_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,wR_arg,wL_arg,Bmag_arg,dBdz_arg):
    tmpwave=((-1)**(m_res_arg+1))*(qi_arg*Ewz_arg*jn((m_res_arg),beta_arg)-wR_arg*pper_arg*jn((m_res_arg+1),beta_arg)+wL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmpadiabat=((pper_arg*pper_arg)/(2*gamma_arg*mi_arg*Bmag_arg))*dBdz_arg
    tmp=tmpwave-tmpadiabat
    return tmp