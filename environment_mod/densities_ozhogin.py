import numpy as np

#####environment_mod.densities_ozhogin############################################

#Description:Routine to calculate the electron density along a magnetic field line after [Ozhogin et al., 2012]
#Inputs:
# L_arg: Lshell
# lamda_arg: geomagnetic latitude in rad
#Outputs:
# neqtmp: equatorial electron number density
# n_lamtmp: electron number density at latitude

#############################################################################

def densities_ozhogin(L_arg,lambda_arg):
    neqtmp=10**(4.4693-0.4903*L_arg)
    lamda_inv=np.arccos(1/L_arg)
    n_lamtmp=neqtmp*((np.cos((np.pi/2)*(lambda_arg/lamda_inv)))**(-0.75))
    return neqtmp,n_lamtmp
