import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def dEkdt(pper_arg,ppar_arg,eta_arg,m_res_arg,Exw_arg,Eyw_arg,Ezw_arg,beta_arg,gamma_arg):
    EwL=0.5*(Exw_arg-Eyw_arg)

    EwR=0.5*(Exw_arg+Eyw_arg)
   
    
    fac=((-1)**(m_res_arg-1))*(const.qe/(gamma_arg*const.me))
    fac1=fac*ppar_arg*Ezw_arg*jn((m_res_arg),beta_arg)
    fac2=-fac*pper_arg*EwL*jn((m_res_arg+1),beta_arg)
    fac3=-fac*pper_arg*EwR*jn((m_res_arg-1),beta_arg)
    
    tmp=(fac1+fac2+fac3)*np.sin(eta_arg)
    return tmp