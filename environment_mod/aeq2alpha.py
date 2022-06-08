import numpy as np
from environment_mod.Bmag_dipole import Bmag_dipole
def aeq2alpha(L_arg,lambda_arg,aeq_arg):
    Blam0=Bmag_dipole(L_arg,lambda_arg)
#    print(Blam0)
    Beq0=Bmag_dipole(L_arg,0)
#    print(Beq0)
    salpha0=np.sin(aeq_arg)*np.sqrt(Blam0/Beq0)
#     print(np.rad2deg(salpha0))
    alpha0=np.arcsin(salpha0)
    
    return alpha0