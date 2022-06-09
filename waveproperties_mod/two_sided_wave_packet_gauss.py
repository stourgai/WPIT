
import numpy as np


def two_sided_wave_packet_gauss(Bw0_arg,lamda_arg,lamda_range):
    Bwave=Bw0_arg*np.exp(-lamda_arg**2/lamda_range**2)

    return Bwave