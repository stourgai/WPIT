import numpy as np
def densities_ozhogin(L_arg,lambda_arg):
    neqtmp=10**(4.4693-0.4903*L_arg)
    lamda_inv=np.arccos(1/L_arg)
    n_lamtmp=neqtmp*((np.cos((np.pi/2)*(lambda_arg/lamda_inv)))**(-0.75))
    return neqtmp,n_lamtmp
