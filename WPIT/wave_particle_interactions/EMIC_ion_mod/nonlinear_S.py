import numpy as np
from scipy.special import jn 

from environment_mod import const

def nonlinear_S(H_arg,wtsq_arg):
    tmp=H_arg/wtsq_arg
    return tmp