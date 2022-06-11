import numpy as np

#####environment_mod.densities_denton############################################

#Description:Routine to calculate the electron density along a magnetic field line after [Denton et al., 2002]
#Inputs:
# ne0_arg: equatorial electron number density
# lamda_arg: geomagnetic latitude in rad
#Outputs:
# n_e_tmp: electron number density

#############################################################################

def densities_denton(ne0_arg,lamda_arg):
    #-----calculate species densities (assuming 94%H+, 5.4%He+, 0.6%O+)
    #lamda geomagnetic latitude
    clat=np.cos(lamda_arg)
    n_e_tmp=ne0_arg*(clat**(-4))

    return n_e_tmp