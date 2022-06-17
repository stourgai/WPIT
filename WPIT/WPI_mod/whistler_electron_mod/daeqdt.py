import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def daeqdt(ppar_arg,pper_arg,alpha_arg,aeq_arg,eta_arg,w1_arg,R1_arg,w2_arg,R2_arg,gamma_arg,beta_arg,wtausq_arg,kz_arg,m_res_arg):
    pmag=np.sqrt(ppar_arg**2+pper_arg**2)
    fac1=-(1/(pmag*pmag))*(np.tan(aeq_arg)/np.tan(alpha_arg))
    fac2a=pow((-1),(m_res_arg-1))*w1_arg*((ppar_arg/gamma_arg)+const.me*R1_arg)*jn((m_res_arg-1),beta_arg)
    fac2b=pow((-1),(m_res_arg-1))*w2_arg*((ppar_arg/gamma_arg)-const.me*R2_arg)*jn((m_res_arg+1),beta_arg)
    fac2=(fac2a-fac2b)*ppar_arg
    fac3=wtausq_arg*const.me*pper_arg/kz_arg
    tmp=fac1*(fac2+fac3)*np.sin(eta_arg)
    
    return tmp
