
import numpy as np

#####waveproperties_mod.two_sided_wave_packet##############################

#Description:Simulate a static, monochromatic and two-sided wave packet
#Inputs:
# Bw0_arg: initial wave B field magnitude
# lamda_arg: latitude in rad
# shape: the higher the sharper the packet edges
# location: logation of HWHM
#Outputs:
# Bwave: wave B field magnitude
##############################################################

def two_sided_wave_packet(Bw0_arg,lamda_arg,shape,location):
    if lamda_arg<0:
        Bwave=Bw0_arg-Bw0_arg*(np.tanh(-1*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    if lamda_arg>0:
        Bwave=Bw0_arg-Bw0_arg*(np.tanh(1*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    return Bwave