
import numpy as np

#####waveproperties_mod.one_sided_wave_packet##############################

#Description:Simulate a static, monochromatic and one-sided wave packet
#Inputs:
# Bw0_arg:Initial wave magnetic field amplitude
# lamda_arg: latitude in rad
# shape: the higher the sharper the packet edges
# location: logation of HWHM
# direction: south or north direction
#Outputs:
# Bwave: wave B field magnitude

#####################################################################

def one_sided_wave_packet(Bw0_arg,lamda_arg,shape,location,direction):
    if direction =='north':
        dir=1
    if direction =='south':
        dir=-1


    Bwave=Bw0_arg*(np.tanh(dir*shape*np.rad2deg(lamda_arg)-2*location)+1)/2
    return Bwave