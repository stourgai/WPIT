import numpy as np
from environment_mod import const
def dwc_ds(wc_arg,lamda_arg,L_arg):

    slat = np.sin(lamda_arg)
    clat = np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    dwce_ds_arg = (3.0*wc_arg/(L_arg*const.Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))
    return dwce_ds_arg