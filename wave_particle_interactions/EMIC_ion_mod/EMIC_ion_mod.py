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

def wpi_params(pper_arg,kper_arg,qi_arg,mi_arg,Bmag_arg,Exw_arg,Eyw_arg,Bxw_arg,Byw_arg,gamma_arg):
    beta_tmp=-(kper_arg*pper_arg)/(qi_arg*Bmag_arg)
    BwR=0.5*(Bxw_arg+Byw_arg)
    BwL=0.5*(Bxw_arg-Byw_arg)
    EwR=0.5*(Exw_arg+Eyw_arg)
    EwL=0.5*(Exw_arg-Eyw_arg)    
    pwR=gamma_arg*mi_arg*(EwR/BwR)
    pwL=gamma_arg*mi_arg*(EwL/BwL)
    wR=(qi_arg*BwR)/(gamma_arg*mi_arg)
    wL=(qi_arg*BwL)/(gamma_arg*mi_arg)    
    
    return beta_tmp,BwR,BwL,EwR,EwL,pwR,pwL,wR,wL

def dzdt(ppar_arg,gamma_arg,mi_arg):
    tmp=ppar_arg/(mi_arg*gamma_arg)
    return tmp

def dppardt(pper_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,wR_arg,wL_arg,Bmag_arg,dBdz_arg):
    tmpwave=((-1)**(m_res_arg+1))*(qi_arg*Ewz_arg*jn((m_res_arg),beta_arg)-wR_arg*pper_arg*jn((m_res_arg+1),beta_arg)+wL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmpadiabat=((pper_arg*pper_arg)/(2*gamma_arg*mi_arg*Bmag_arg))*dBdz_arg
    tmp=tmpwave-tmpadiabat
    return tmp

def dpperdt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg,Bmag_arg,dBdz_arg):
    tmpwave=((-1)**(m_res_arg+1))*((ppar_arg-pwR_arg)*wR_arg*jn((m_res_arg+1),beta_arg)-(ppar_arg-pwL_arg)*wL_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmpadiabat=((pper_arg*ppar_arg)/(2*gamma_arg*mi_arg*Bmag_arg))*dBdz_arg
    tmp=tmpwave+tmpadiabat
    return tmp

def detadt(ppar_arg,mres_arg,wc_arg,gamma_arg,kpar_arg,mi_arg,w_arg):
    tmp=((mres_arg*wc_arg)/gamma_arg)+((kpar_arg*ppar_arg)/(gamma_arg*mi_arg))-w_arg
    return tmp

def dlamdadt(ppar_arg,lamda_arg,gamma_arg,mi_arg,L_arg):
    ork=ppar_arg/(gamma_arg*mi_arg*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork

def dalphadt(pper_arg,ppar_arg,eta_arg,Ewz_arg,m_res_arg,qi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg):
    fac=((-1)**(m_res_arg+1))
    pmagsq=ppar_arg**2+pper_arg**2
    fac1=(fac/pmagsq)*np.sin(eta_arg)
    fac2a=ppar_arg*qi_arg*Ewz_arg*jn((m_res_arg),beta_arg)
    fac2b=pper_arg*pwR_arg*wR_arg*jn((m_res_arg+1),beta_arg)
    fac2c=pper_arg*pwL_arg*wL_arg*jn((m_res_arg-1),beta_arg)
    tmp=fac1*(fac2a-fac2b-fac2c)
    return tmp

def daeqdt(pper_arg,ppar_arg,eta_arg,aeq_arg,Ewz_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,pwR_arg,pwL_arg,beta_arg,wR_arg,wL_arg,Bwytmp_arg,wmega_arg,kappa_arg):
    fac1=((-1)**(m_res_arg+1))
    pmagsq=ppar_arg**2+pper_arg**2
    fac2=(ppar_arg/pper_arg)*((np.tan(aeq_arg)*np.sin(eta_arg))/pmagsq)
    fac3=-qi_arg*Ewz_arg*pper_arg*jn((m_res_arg),beta_arg)+wR_arg*(pmagsq-pwR_arg*ppar_arg)*jn((m_res_arg+1),beta_arg)-wL_arg*(pmagsq+pwL_arg*ppar_arg)*jn((m_res_arg-1),beta_arg)
    tmp=fac1*fac2*fac3
    return tmp

def dgammadt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,EwR_arg,EwL_arg):
    fac1=((-1)**(m_res_arg+1))*(qi_arg/(mi_arg*mi_arg*const.c_light*const.c_light*gamma_arg))
    fac2=(Ewz_arg*ppar_arg*jn((m_res_arg),beta_arg)-EwR_arg*pper_arg*jn((m_res_arg+1),beta_arg)-EwL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp

def dEkdt(pper_arg,ppar_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,EwR_arg,EwL_arg,wmega_arg,kappa_arg):
    fac1=((-1)**(m_res_arg+1))*(qi_arg/(mi_arg*gamma_arg))
    fac2=(Ewz_arg*ppar_arg*jn((m_res_arg),beta_arg)-EwR_arg*pper_arg*jn((m_res_arg+1),beta_arg)-EwL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmp=fac1*fac2
    return tmp

def nonlinear_S(H_arg,wtsq_arg):
    tmp=H_arg/wtsq_arg
    return tmp
    
def nonlinear_H(pper_arg,ppar_arg,kpar_arg,gamma_arg,mres_arg,mi_arg,wce_arg,dkpar_dt_arg,dwcdz_arg,dwdt_arg):
    dwc_dt=(ppar_arg/(gamma_arg*mi_arg))*dwcdz_arg
    fac1=(mres_arg/gamma_arg)*dwc_dt
    fac2=(ppar_arg/(gamma_arg*mi_arg))*dkpar_dt_arg
    fac3=-((kpar_arg*pper_arg*pper_arg)/(2*gamma_arg*gamma_arg*mi_arg*mi_arg*wce_arg))*dwcdz_arg
    tmp=fac1+fac2+fac3-dwdt_arg
    return tmp

def nonlinear_theta(C0_arg,Cp1_arg,Cm1_arg,m_res_arg,beta_arg):
    tmp=C0_arg*jn((m_res_arg),beta_arg)+Cp1_arg*jn((m_res_arg+1),beta_arg)+Cm1_arg*jn((m_res_arg-1),beta_arg)
    tmpwtsq=np.abs(tmp)
    return tmp,tmpwtsq

def nonlinear_C0(ppar_arg,kpar_arg,mres_arg,gamma_arg,qi_arg,mi_arg,wce_arg,Ewz_arg):
    fac1=(-1)**(mres_arg+1)
    fac2a=(qi_arg*kpar_arg)/(gamma_arg*mi_arg)
    fac2b=((qi_arg*ppar_arg)/(gamma_arg*gamma_arg*gamma_arg*mi_arg*mi_arg*const.c_light*const.c_light))*(mres_arg*wce_arg+((kpar_arg*ppar_arg)/mi_arg))
    tmp=fac1*(fac2a-fac2b)*Ewz_arg
    return tmp

def nonlinear_C1p(pper_arg,ppar_arg,kpar_arg,mres_arg,qi_arg,mi_arg,gamma_arg,wR_arg,EwR_arg,wce_arg):
    fac1=(-1)**(mres_arg+1)
    fac2a=-(wR_arg*kpar_arg)/(gamma_arg*mi_arg)
    fac2b=((qi_arg*EwR_arg)/(gamma_arg*gamma_arg*gamma_arg*mi_arg*mi_arg*const.c_light*const.c_light))*(mres_arg*wce_arg+((kpar_arg*ppar_arg)/mi_arg))
    tmp=fac1*(fac2a+fac2b)*pper_arg
    return tmp

def nonlinear_C1m(pper_arg,ppar_arg,kpar_arg,mres_arg,qi_arg,mi_arg,gamma_arg,wL_arg,EwL_arg,wce_arg):
    fac1=(-1)**(mres_arg+1)
    fac2a=(wL_arg*kpar_arg)/(gamma_arg*mi_arg)
    fac2b=((qi_arg*EwL_arg)/(gamma_arg*gamma_arg*gamma_arg*mi_arg*mi_arg*const.c_light*const.c_light))*(mres_arg*wce_arg+((kpar_arg*ppar_arg)/mi_arg))
    tmp=fac1*(fac2a+fac2b)*pper_arg
    return tmp

def dwcdt(ppar,m,gamma,dwcdz):
    tmp=(ppar/(m*gamma))*dwcdz
    return tmp

def dkpardt(ppar,m,gamma,psi,dkdz):
    tmp=(ppar/(m*gamma))*dkdz*np.cos(psi)
    return tmp