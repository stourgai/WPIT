import numpy as np
def omega_lh(wce_arg,wpe_arg,wci_arg,wpi_arg):
    fac1=1/(wci_arg*wce_arg)
    fac2=1/(wpe_arg**2+wpi_arg**2)
    sqrtfac=np.sqrt(fac1+fac2)
    wlh_tmp=1/sqrtfac
    return wlh_tmp