import numpy as np

#####environment_mod.omega_uh############################################

#Description:Routine to calculate the upper hybrid resonance frequency
#Inputs:
# wce_arg: electron gyrofrequency
# wpe_arg: electron plasma frequency 
#Outputs:
# tmp: upper hybrid resonance frequency

#############################################################################

def omega_uh(wce_arg,wpe_arg):
    tmpsq=wpe_arg*wpe_arg+wce_arg*wce_arg
    tmp=np.sqrt(tmpsq)
    return tmp