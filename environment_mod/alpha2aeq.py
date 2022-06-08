import numpy as np
from environment_mod.Bmag_dipole import Bmag_dipole

def alpha2aeq(L_arg,lambda_arg,alpha_arg):
    Blam0=Bmag_dipole(L_arg,lambda_arg)
    Beq0=Bmag_dipole(L_arg,0)
    salphaeq0=np.sin(alpha_arg)*np.sqrt(Beq0/Blam0)

    alphaeq0=np.arcsin(salphaeq0)
    return alphaeq0