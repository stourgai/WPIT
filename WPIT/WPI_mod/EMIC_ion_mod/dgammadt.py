import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def dgammadt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,EwR_arg,EwL_arg):
    fac1=((-1)**(m_res_arg+1))*(qi_arg/(mi_arg*mi_arg*const.c_light*const.c_light*gamma_arg))
    fac2=(Ewz_arg*ppar_arg*jn((m_res_arg),beta_arg)-EwR_arg*pper_arg*jn((m_res_arg+1),beta_arg)-EwL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp