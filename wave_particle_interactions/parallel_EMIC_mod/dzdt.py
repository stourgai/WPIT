import numpy as np
from scipy.special import jn 

from environment_mod import const

def dzdt(ppar_arg,gamma_arg,mi_arg):
    krk=ppar_arg/(gamma_arg*mi_arg)
    return krk
