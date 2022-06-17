import numpy as np
from scipy.special import jn 

from environment_mod import const

def nonlinear_S(H,wtsq):
    tmp=H/wtsq
    return tmp