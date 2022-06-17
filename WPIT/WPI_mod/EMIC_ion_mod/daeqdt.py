import numpy as np
from scipy.special import jn 

def daeqdt(pper_arg,ppar_arg,eta_arg,aeq_arg,Ewz_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg):
    fac1=((-1)**(m_res_arg+1))
    pmagsq=ppar_arg**2+pper_arg**2
    fac2=(ppar_arg/pper_arg)*((np.tan(aeq_arg)*np.sin(eta_arg))/pmagsq)
    fac3=-qi_arg*Ewz_arg*pper_arg*jn((m_res_arg),beta_arg)+wR_arg*(pmagsq-pwR_arg*ppar_arg)*jn((m_res_arg+1),beta_arg)-wL_arg*(pmagsq+pwL_arg*ppar_arg)*jn((m_res_arg-1),beta_arg)
    tmp=fac1*fac2*fac3
    return tmp