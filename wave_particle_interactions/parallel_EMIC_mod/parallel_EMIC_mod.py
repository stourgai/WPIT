import numpy as np
from scipy.special import jn 

import os
import sys
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "../../../../environment")
sys.path.append(fpath)
# print(current_dir)
# print(fpath)
import  const

def dzdt(ppar_arg,gamma_arg,mi_arg):
    krk=ppar_arg/(gamma_arg*mi_arg)
    return krk

def dppardt(pper_arg,eta_arg,wc_arg,Bw_arg,gamma_arg,dwcds_arg,q_arg,m_arg):
    lrk1=((q_arg*Bw_arg)/(gamma_arg*m_arg))*pper_arg*np.sin(eta_arg)
    lrk2=((pper_arg*pper_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    lrk=lrk1-lrk2
    return lrk

def dpperdt(ppar_arg,pper_arg,eta_arg,Bw_arg,gamma_arg,wmega_arg,kappa_arg,wc_arg,dwcds_arg,q_arg,m_arg):
    mrk1=(q_arg*Bw_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.sin(eta_arg)
    mrk2=((pper_arg*ppar_arg)/(2*m_arg*gamma_arg*wc_arg))*dwcds_arg
    mrk=mrk1+mrk2
    return mrk

def detadt(ppar_arg,pper_arg,eta_arg,Bw_arg,wmega_arg,kappa_arg,wc_arg,gamma_arg,q_arg,m_arg):
    nrk1=((q_arg*Bw_arg)/pper_arg)*((wmega_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*np.cos(eta_arg)
    nrk2a=(kappa_arg*ppar_arg)/(gamma_arg*m_arg)
    nrk2b=-wmega_arg
    nrk2c=-wc_arg/gamma_arg
    nrk=nrk1+nrk2a+nrk2b+nrk2c
    return nrk

def dlamdadt(ppar_arg,lamda_arg,gamma_arg,m_arg,L_arg):
    ork=ppar_arg/(gamma_arg*m_arg*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork

def dalphadt(ppar_arg,pper_arg,eta_arg,gamma_arg,wce_arg,dwcds_arg,w_wave_arg,kappa_arg,Bw_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    a_rk1=((q_arg*Bw_arg)/(p_mag**2))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-((pper_arg**2)/(gamma_arg*m_arg)))*np.sin(eta_arg)
    a_rk2=(pper_arg/(2*gamma_arg*m_arg*wce_arg))*dwcds_arg
    a_rk=a_rk1+a_rk2
    return a_rk

def daeqdt(ppar_arg,pper_arg,eta_arg,gamma_arg,w_wave_arg,kappa_arg,Bw_arg,aeq_arg,alpha_arg,m_arg,q_arg):
    p_mag=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
    aeq_rk=((q_arg*Bw_arg)/(p_mag**2))*(np.tan(aeq_arg)/np.tan(alpha_arg))*(((w_wave_arg/kappa_arg)-(ppar_arg/(gamma_arg*m_arg)))*ppar_arg-(pper_arg**2/(gamma_arg*m_arg)))*np.sin(eta_arg)
    return aeq_rk

def dgammadt(pper_arg,eta_arg,gamma_arg,Bw_arg,kappa_arg,w_arg,q_arg,m_arg):
    fac1=1/(gamma_arg*m_arg*m_arg*const.c_light*const.c_light)
    fac2=q_arg*pper_arg*Bw_arg*(w_arg/kappa_arg)*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp

def dEkdt(pper_arg,eta_arg,gamma_arg,Bw_arg,kappa_arg,w_arg,q_arg,m_arg):
    fac1=1/(gamma_arg*m_arg)
    fac2=q_arg*pper_arg*Bw_arg*(w_arg/kappa_arg)*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp

def nonlinear_S(H_arg,wtsq_arg):
    tmp=H_arg/wtsq_arg
    return tmp

def nonlinear_H(pper_arg,ppar_arg,kappa_arg,gamma_arg,m_arg,wce_arg,dk_dt_arg,dwcdz_arg,dwdt_arg):
    dwc_dt=(ppar_arg/(gamma_arg*m_arg))*dwcdz_arg
    fac1=-(1/gamma_arg)*dwc_dt
    fac2=(ppar_arg/(gamma_arg*m_arg))*dk_dt_arg
    fac3=-((kappa_arg*pper_arg*pper_arg)/(2*gamma_arg*gamma_arg*m_arg*m_arg*wce_arg))*dwcdz_arg
    tmp=fac1+fac2+fac3-dwdt_arg
    return tmp

def nonlinear_theta(pper_arg,ppar_arg,Bw_arg,kappa_arg,gamma_arg,m_arg,q_arg,wce_arg,w_arg):
    fac1=q_arg*Bw_arg*pper_arg
    fac2a=kappa_arg/(gamma_arg*gamma_arg*m_arg*m_arg)
    fac2b=(wce_arg-((kappa_arg*ppar_arg)/m_arg))*(w_arg/(kappa_arg*gamma_arg*gamma_arg*gamma_arg*m_arg*const.c_light*const.c_light))
    tmp=fac1*(fac2a+fac2b)
    tmpwtsq_arg=np.abs(tmp)
    return tmp,tmpwtsq_arg