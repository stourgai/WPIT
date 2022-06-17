import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def dalphadt(pper_arg,ppar_arg,eta_arg,w1_arg,w2_arg,R1_arg,R2_arg,wtau_sq_arg,kz_arg,beta_arg,m_res_arg,gamma_arg,wce_arg,dwhds_arg):
    pmag=np.sqrt(pper_arg*pper_arg+ppar_arg*ppar_arg)
    fac1=(1/(pmag*pmag))
    mrk=(-pow((-1),(m_res_arg-1))*(w1_arg*((ppar_arg/gamma_arg)+(const.me*R1_arg))*jn((m_res_arg-1),beta_arg)-
            w2_arg*((ppar_arg/gamma_arg)-(const.me*R2_arg))*jn((m_res_arg+1),beta_arg))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg*ppar_arg)/(2*wce_arg))*dwhds_arg)
    fac2a=ppar_arg*mrk
    
    lrk=((wtau_sq_arg*const.me)/kz_arg)*np.sin(eta_arg)-(1/(const.me*gamma_arg))*((pper_arg*pper_arg)/(2*wce_arg))*dwhds_arg
    fac2b=pper_arg*lrk
    
    tmp=fac1*(fac2a-fac2b)
    
    return tmp
