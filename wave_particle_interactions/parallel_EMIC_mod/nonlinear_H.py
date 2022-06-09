import numpy as np
from scipy.special import jn 

from environment_mod import const

def nonlinear_H(pper_arg,ppar_arg,kappa_arg,gamma_arg,m_arg,wce_arg,dk_dt_arg,dwcdz_arg,dwdt_arg):
    dwc_dt=(ppar_arg/(gamma_arg*m_arg))*dwcdz_arg
    fac1=-(1/gamma_arg)*dwc_dt
    fac2=(ppar_arg/(gamma_arg*m_arg))*dk_dt_arg
    fac3=-((kappa_arg*pper_arg*pper_arg)/(2*gamma_arg*gamma_arg*m_arg*m_arg*wce_arg))*dwcdz_arg
    tmp=fac1+fac2+fac3-dwdt_arg
    return tmp