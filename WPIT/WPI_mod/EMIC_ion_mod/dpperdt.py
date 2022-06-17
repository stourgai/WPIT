import numpy as np
from scipy.special import jn 


def dpperdt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg,Bmag_arg,dBdz_arg):
    tmpwave=((-1)**(m_res_arg+1))*((ppar_arg-pwR_arg)*wR_arg*jn((m_res_arg+1),beta_arg)-(ppar_arg-pwL_arg)*wL_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmpadiabat=((pper_arg*ppar_arg)/(2*gamma_arg*mi_arg*Bmag_arg))*dBdz_arg
    tmp=tmpwave+tmpadiabat
    return tmp