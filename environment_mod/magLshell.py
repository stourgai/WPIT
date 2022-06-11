import numpy as np
from environment_mod import const

#####environment_mod.magLshell############################################

#Description:Calculate the Lshell given the latitude and the geocentric distance
#Inputs:
# r_arg: geocentric distance in km
# lat_arg: magnetic latitude in rad
#Outputs:
# L_tmp: L shell

#############################################################################

def magLshell(r_arg,lat_arg):
    L_tmp=r_arg/(const.Re*np.cos(lat_arg)*np.cos(lat_arg))
    return L_tmp