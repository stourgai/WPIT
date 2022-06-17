import numpy as np


def nonlinear_H(pper,ppar,kpar,gamma,m_res,m,wce,dkpar_dt,dwcdz,dwdt):
    dwc_dt=(ppar/(gamma*m))*dwcdz
    fac1=(m_res/gamma)*dwc_dt
    fac2=-(ppar/(gamma*m))*dkpar_dt
    fac3=((kpar*pper*pper)/(2*gamma*gamma*m*m*wce))*dwcdz
    tmp=fac1+fac2+fac3-dwdt
    return tmp