import numpy as np
from scipy.special import jn 

def dEkdt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,EwR_arg,EwL_arg,wmega_arg,kappa_arg):
    fac1=((-1)**(m_res_arg+1))*(qi_arg/(mi_arg*gamma_arg))
    fac2=(Ewz_arg*ppar_arg*jn((m_res_arg),beta_arg)-EwR_arg*pper_arg*jn((m_res_arg+1),beta_arg)-EwL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp