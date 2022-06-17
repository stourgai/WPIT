import numpy as np

from WPIT.Environment_mod import const

def nonlinear_C0(ppar_arg,kpar_arg,mres_arg,gamma_arg,qi_arg,mi_arg,wce_arg,Ewz_arg):
    fac1=(-1)**(mres_arg+1)
    fac2a=(qi_arg*kpar_arg)/(gamma_arg*mi_arg)
    fac2b=((qi_arg*ppar_arg)/(gamma_arg*gamma_arg*gamma_arg*mi_arg*mi_arg*const.c_light*const.c_light))*(mres_arg*wce_arg+((kpar_arg*ppar_arg)/mi_arg))
    tmp=fac1*(fac2a-fac2b)*Ewz_arg
    return tmp