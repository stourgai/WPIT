import numpy as np
from scipy.special import jn 

from environment_mod import const

def wpi_params(pper_arg,kper_arg,qi_arg,mi_arg,Bmag_arg,Exw_arg,Eyw_arg,Bxw_arg,Byw_arg,gamma_arg):
    beta_tmp=-(kper_arg*pper_arg)/(qi_arg*Bmag_arg)
    BwR=0.5*(Bxw_arg+Byw_arg)
    BwL=0.5*(Bxw_arg-Byw_arg)
    EwR=0.5*(Exw_arg+Eyw_arg)
    EwL=0.5*(Exw_arg-Eyw_arg)    
    pwR=gamma_arg*mi_arg*(EwR/BwR)
    pwL=gamma_arg*mi_arg*(EwL/BwL)
    wR=(qi_arg*BwR)/(gamma_arg*mi_arg)
    wL=(qi_arg*BwL)/(gamma_arg*mi_arg)    
    
    return beta_tmp,BwR,BwL,EwR,EwL,pwR,pwL,wR,wL
