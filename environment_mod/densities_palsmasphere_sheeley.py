import numpy as np

#####environment_mod.densities_palsmasphere_sheeley############################################

#Description:Routine to calculate the gradient of magnetic field strength with respect to the distance along the magnetic field line
#Inputs:
# B_arg: magnetic field in T
# lamda_arg: magnetic latitude in rad
# L_arg: L shell
#Outputs:
# dB_ds_arg:  gradient of magnetic field along the field line


#############################################################################

def densities_palsmasphere_sheeley(L_arg):
    ne_mean=1390*(3/L_arg)**4.8
    ne_min=ne_mean-440*(3/L_arg)**3.6
    ne_max=ne_mean+440*(3/L_arg)**3.6

    return ne_mean,ne_min,ne_max
