import numpy as np
from environment_mod import const
def debye_length(ne_arg,Te_arg):
    debye_fac_1=const.epsilon0*const.kb*Te_arg
    debye_fac_2=ne_arg*const.qe*const.qe
    tmp=debye_fac_1/debye_fac_2

    return tmp
    