import numpy as np
from scipy.special import jn 

def dkpardt(ppar,m,gamma,psi,dkdz):
    tmp=(ppar/(m*gamma))*dkdz*np.cos(psi)
    return tmp