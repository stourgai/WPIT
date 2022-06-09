import numpy as np
from scipy.special import jn 

from environment_mod import const


def dzdt(ppar_arg,gamma_arg,mi_arg):
    tmp=ppar_arg/(mi_arg*gamma_arg)
    return tmp