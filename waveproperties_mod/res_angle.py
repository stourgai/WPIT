import numpy as np

#####waveproperties_mod.res_angle##############################

#Description:Routine to calculate the resonance angle
#Inputs:
# P_arg: Stix P parameter
# S_arg: Stix S parameter
#Outputs:
# thetares: resonance angle in rad

##############################################################

def res_angle(P_arg,S_arg):
    tansq=-P_arg/S_arg
    tan=np.sqrt(tansq)
    thetares=np.arctan(tan)
    return thetares

