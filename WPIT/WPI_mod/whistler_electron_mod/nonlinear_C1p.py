import numpy as np


from WPIT.Environment_mod import const

def nonlinear_C1p(pper_arg,ppar_arg,w2_arg,Exw_arg,Eyw_arg,m_res_arg,wce_arg,kz_arg,gamma_arg):
    EwL=0.5*(Exw_arg-Eyw_arg)
    tau=m_res_arg*wce_arg-(kz_arg*ppar_arg/const.me)
    fac_mres=(-1)**(m_res_arg-1)

    fac1a=(tau*const.qe*pper_arg*EwL)/(gamma_arg*gamma_arg*gamma_arg*const.me*const.me*const.c_light*const.c_light)
    fac1b=(pper_arg*kz_arg*w2_arg)/(gamma_arg*gamma_arg*const.me)

    tmp=fac_mres*(fac1a+fac1b)
    return tmp