import numpy as np
from environment_mod import const
def magLshell(r_arg,lat_arg):
    L_tmp=r_arg/(const.Re*np.cos(lat_arg)*np.cos(lat_arg))
    return L_tmp