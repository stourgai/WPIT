import numpy as np
from scipy.special import jn 

def dalphadt(pper_arg,ppar_arg,eta_arg,Ewz_arg,m_res_arg,qi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg):
    fac=((-1)**(m_res_arg+1))
    pmagsq=ppar_arg**2+pper_arg**2
    fac1=(fac/pmagsq)*np.sin(eta_arg)
    fac2a=ppar_arg*qi_arg*Ewz_arg*jn((m_res_arg),beta_arg)
    fac2b=pper_arg*pwR_arg*wR_arg*jn((m_res_arg+1),beta_arg)
    fac2c=pper_arg*pwL_arg*wL_arg*jn((m_res_arg-1),beta_arg)
    tmp=fac1*(fac2a-fac2b-fac2c)
    return tmp