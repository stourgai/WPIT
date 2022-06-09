import numpy as np
from scipy.special import jn 

from environment_mod import const

def detadt(ppar_arg,mres_arg,wc_arg,gamma_arg,kpar_arg,mi_arg,w_arg):
    tmp=((mres_arg*wc_arg)/gamma_arg)+((kpar_arg*ppar_arg)/(gamma_arg*mi_arg))-w_arg
    return tmp