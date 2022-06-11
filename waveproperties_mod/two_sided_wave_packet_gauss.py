
import numpy as np

#####waveproperties_mod.two_sided_wave_packet_gauss##############

#Description:Simulate a static, monochromatic gaussian wave packet
#Inputs:
# Bw0_arg: initial wave B field magnitude
# lamda_arg: latitude in rad
# lamdaw_arg: latitudinal range
#Outputs:
# Bwave: wave B field magnitude
#################################################################

def two_sided_wave_packet_gauss(Bw0_arg,lamda_arg,lamda_range):
    Bwave=Bw0_arg*np.exp(-lamda_arg**2/lamda_range**2)

    return Bwave