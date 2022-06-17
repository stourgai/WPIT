import numpy as np
from WPIT.Environment_mod import const

def nonlinear_C0(ppar,m_res,wce,kz,gamma,Ezw):
    tau=m_res*wce-(kz*ppar/const.me)
    facres=(-1)**(m_res-1)
    fac1=const.qe*kz/(gamma*const.me)
    fac2=(tau*const.qe*ppar)/(gamma*gamma*gamma*const.me*const.me*const.c_light*const.c_light)
    tmp=-facres*(fac1+fac2)*Ezw
    
    return tmp
