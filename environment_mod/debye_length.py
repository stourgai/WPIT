import numpy as np
from environment_mod import const

#####environment_mod.debye_length############################################

#Description:Routine to calculate the debye length
#Inputs:
# ne_arg: electron density in m^-3
# Te_arg: electron temperature in K
#Outputs:
# tmp: Debye length in m

#############################################################################

def debye_length(ne_arg,Te_arg):
    debye_fac_1=const.epsilon0*const.kb*Te_arg
    debye_fac_2=ne_arg*const.qe*const.qe
    tmp=debye_fac_1/debye_fac_2

    return tmp
    