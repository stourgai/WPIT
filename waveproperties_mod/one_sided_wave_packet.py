
import numpy as np


def one_sided_wave_packet(Bw0_arg,lamda_arg,shape,location,direction):
    if direction =='north':
        dir=1
    if direction =='south':
        dir=-1


    Bwave=Bw0_arg*(np.tanh(dir*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    return Bwave