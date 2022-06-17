import numpy as np

from WPIT.Environment_mod import const

def nonlinear_C1p(pper_arg,ppar_arg,kpar_arg,mres_arg,qi_arg,mi_arg,gamma_arg,wR_arg,EwR_arg,wce_arg):
    fac1=(-1)**(mres_arg+1)
    fac2a=-(wR_arg*kpar_arg)/(gamma_arg*mi_arg)
    fac2b=((qi_arg*EwR_arg)/(gamma_arg*gamma_arg*gamma_arg*mi_arg*mi_arg*const.c_light*const.c_light))*(mres_arg*wce_arg+((kpar_arg*ppar_arg)/mi_arg))
    tmp=fac1*(fac2a+fac2b)*pper_arg
    return tmp