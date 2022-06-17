import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

def detadt(ppar_arg,m_res_arg,wce_arg,wwave_arg,gamma_arg,kz_arg):
    nrk=((m_res_arg*wce_arg)/gamma_arg)-wwave_arg-kz_arg*(ppar_arg/(const.me*gamma_arg))
    return nrk