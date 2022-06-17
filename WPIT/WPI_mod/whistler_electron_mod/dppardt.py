import numpy as np

from WPIT.Environment_mod import const

def dppardt(pper_arg,eta_arg,wtau_sq_arg,kz_arg,gamma_arg,wce_arg,dwds_arg):
    lrk=((wtau_sq_arg*const.me)/kz_arg)*np.sin(eta_arg)-(1/(const.me*gamma_arg))*((pper_arg*pper_arg)/(2*wce_arg))*dwds_arg
    return lrk     