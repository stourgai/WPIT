import numpy as np

from WPIT.Environment_mod import const

def dlamdadt(ppar_arg,lamda_arg,gamma_arg,mi_arg,L_arg):
    ork=ppar_arg/(gamma_arg*mi_arg*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork