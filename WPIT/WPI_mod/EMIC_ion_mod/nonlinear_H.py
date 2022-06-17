import numpy as np

def nonlinear_H(pper_arg,ppar_arg,kpar_arg,gamma_arg,mres_arg,mi_arg,wce_arg,dkpar_dt_arg,dwcdz_arg,dwdt_arg):
    dwc_dt=(ppar_arg/(gamma_arg*mi_arg))*dwcdz_arg
    fac1=(mres_arg/gamma_arg)*dwc_dt
    fac2=(ppar_arg/(gamma_arg*mi_arg))*dkpar_dt_arg
    fac3=-((kpar_arg*pper_arg*pper_arg)/(2*gamma_arg*gamma_arg*mi_arg*mi_arg*wce_arg))*dwcdz_arg
    tmp=fac1+fac2+fac3-dwdt_arg
    return tmp