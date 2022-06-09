import numpy as np
from scipy.special import jn 

from environment_mod import const

def dwcdt(ppar,m,gamma,dwcdz):
    tmp=(ppar/(m*gamma))*dwcdz
    return tmp