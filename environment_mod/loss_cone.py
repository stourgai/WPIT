import numpy as np
def loss_cone(L_arg):
    fac=1/np.sqrt(L_arg**3*np.sqrt(4-3/L_arg))
    tmp=np.arcsin(fac)

    return tmp
