import numpy as np
def omega_uh(wce_arg,wpe_arg):
    tmpsq=wpe_arg*wpe_arg+wce_arg*wce_arg
    tmp=np.sqrt(tmpsq)
    return tmp