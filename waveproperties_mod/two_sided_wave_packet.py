
import numpy as np


def two_sided_wave_packet(Bw0_arg,lamda_arg,shape,location):
    if lamda_arg<0:
        Bwave=Bw0_arg-Bw0_arg*(np.tanh(-1*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    if lamda_arg>0:
        Bwave=Bw0_arg-Bw0_arg*(np.tanh(1*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    return Bwave