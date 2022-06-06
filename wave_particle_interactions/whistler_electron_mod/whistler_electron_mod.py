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

def wpi_params(m_res_arg,ppar_arg,pper_arg,Bxw_arg, Byw_arg,Exw_arg,Eyw_arg,Ezw_arg,kz_arg,kx_arg,wce_arg):
    p_mag_arg=np.sqrt(ppar_arg*ppar_arg+pper_arg*pper_arg)
       
    gamma_arg=np.sqrt(1+(p_mag_arg/(const.me*const.c_light))*(p_mag_arg/(const.me*const.c_light)))
    w1_arg=(const.qe/(2*const.me))*(Bxw_arg+Byw_arg)   #Borntik thesis 2.25e
    w2_arg=(const.qe/(2*const.me))*(Bxw_arg-Byw_arg)   #Borntik thesis 2.25e
    wtau0_sq_arg=(w1_arg*kz_arg*pper_arg)/(gamma_arg*const.me)   #Borntik thesis 2.25d
    beta_arg=(kx_arg*pper_arg)/(const.me*gamma_arg*wce_arg)   #Borntik thesis 2.25a
    a1_arg=w2_arg/w1_arg   #Borntik thesis 2.25f
    a2_arg=(const.qe*Ezw_arg)/(w1_arg*pper_arg)   #Borntik thesis 2.25g
    R1_arg=(Exw_arg+Eyw_arg)/(Bxw_arg+Byw_arg)   #Borntik thesis 2.25h
    R2_arg=(Exw_arg-Eyw_arg)/(Bxw_arg-Byw_arg)   #Borntik thesis 2.25h

    wtau_sq_arg = (pow((-1),(m_res_arg-1)) * wtau0_sq_arg * 
            ( jn( (m_res_arg-1), beta_arg ) - 
                a1_arg*jn( (m_res_arg+1) , beta_arg ) +
                gamma_arg*a2_arg*jn( m_res_arg , beta_arg ) ))    #Borntik thesis 2.25c

    return(gamma_arg,w1_arg,w2_arg,wtau_sq_arg,R1_arg,R2_arg,beta_arg)


def dzdt(ppar_arg,gamma_arg):
    krk=ppar_arg/(gamma_arg*const.me)
    return krk

def dppardt(pper_arg,eta_arg,wtau_sq_arg,kz_arg,gamma_arg,wce_arg,dwds_arg):
    lrk=((wtau_sq_arg*const.me)/kz_arg)*np.sin(eta_arg)-(1/(const.me*gamma_arg))*((pper_arg*pper_arg)/(2*wce_arg))*dwds_arg
    return lrk    

def dpperdt(ppar_arg,pper_arg,eta_arg,w1_arg,w2_arg,beta_arg,gamma_arg,R1_arg,R2_arg,m_res_arg,wce_arg,dwds_arg):
    mrk=(-pow((-1),(m_res_arg-1))*(w1_arg*((ppar_arg/gamma_arg)+(const.me*R1_arg))*jn((m_res_arg-1),beta_arg)-
            w2_arg*((ppar_arg/gamma_arg)-(const.me*R2_arg))*jn((m_res_arg+1),beta_arg))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg*ppar_arg)/(2*wce_arg))*dwds_arg)
    return mrk

def detadt(ppar_arg,m_res_arg,wce_arg,wwave_arg,gamma_arg,kz_arg):
    nrk=((m_res_arg*wce_arg)/gamma_arg)-wwave_arg-kz_arg*(ppar_arg/(const.me*gamma_arg))
    return nrk

def dlamdadt(ppar_arg,lamda_arg,gamma_arg,L_arg):
    ork=ppar_arg/(gamma_arg*const.me*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork

def dalphadt(pper_arg,ppar_arg,eta_arg,w1_arg,w2_arg,R1_arg,R2_arg,wtau_sq_arg,kz_arg,beta_arg,m_res_arg,gamma_arg,wce_arg,dwhds_arg):
    pmag=np.sqrt(pper_arg*pper_arg+ppar_arg*ppar_arg)
    fac1=(1/(pmag*pmag))
    mrk=(-pow((-1),(m_res_arg-1))*(w1_arg*((ppar_arg/gamma_arg)+(const.me*R1_arg))*jn((m_res_arg-1),beta_arg)-
            w2_arg*((ppar_arg/gamma_arg)-(const.me*R2_arg))*jn((m_res_arg+1),beta_arg))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg*ppar_arg)/(2*wce_arg))*dwhds_arg)
    fac2a=ppar_arg*mrk
    
    lrk=((wtau_sq_arg*const.me)/kz_arg)*np.sin(eta_arg)-(1/(const.me*gamma_arg))*((pper_arg*pper_arg)/(2*wce_arg))*dwhds_arg
    fac2b=pper_arg*lrk
    
    tmp=fac1*(fac2a-fac2b)
    
    return tmp


def daeqdt(ppar_arg,pper_arg,alpha_arg,aeq_arg,eta_arg,w1_arg,R1_arg,w2_arg,R2_arg,gamma_arg,beta_arg,wtausq_arg,kz_arg,m_res_arg):
    pmag=np.sqrt(ppar_arg**2+pper_arg**2)
    fac1=-(1/(pmag*pmag))*(np.tan(aeq_arg)/np.tan(alpha_arg))
    fac2a=pow((-1),(m_res_arg-1))*w1_arg*((ppar_arg/gamma_arg)+const.me*R1_arg)*jn((m_res_arg-1),beta_arg)
    fac2b=pow((-1),(m_res_arg-1))*w2_arg*((ppar_arg/gamma_arg)-const.me*R2_arg)*jn((m_res_arg+1),beta_arg)
    fac2=(fac2a-fac2b)*ppar_arg
    fac3=wtausq_arg*const.me*pper_arg/kz_arg
    tmp=fac1*(fac2+fac3)*np.sin(eta_arg)
    
    return tmp

def dgammadt(pper_arg,ppar_arg,eta_arg,m_res_arg,Exw_arg,Eyw_arg,Ezw_arg,beta_arg,gamma_arg):
    EwL=0.5*(Exw_arg-Eyw_arg)

    EwR=0.5*(Exw_arg+Eyw_arg)
   
    
    fac=((-1)**(m_res_arg-1))*(const.qe/(gamma_arg*const.me*const.me*const.c_light*const.c_light))
    fac1=fac*ppar_arg*Ezw_arg*jn((m_res_arg),beta_arg)
    fac2=-fac*pper_arg*EwL*jn((m_res_arg+1),beta_arg)
    fac3=-fac*pper_arg*EwR*jn((m_res_arg-1),beta_arg)
    
    tmp=(fac1+fac2+fac3)*np.sin(eta_arg)
    return tmp

def dEkdt(pper_arg,ppar_arg,eta_arg,m_res_arg,Exw_arg,Eyw_arg,Ezw_arg,beta_arg,gamma_arg):
    EwL=0.5*(Exw_arg-Eyw_arg)

    EwR=0.5*(Exw_arg+Eyw_arg)
   
    
    fac=((-1)**(m_res_arg-1))*(const.qe/(gamma_arg*const.me))
    fac1=fac*ppar_arg*Ezw_arg*jn((m_res_arg),beta_arg)
    fac2=-fac*pper_arg*EwL*jn((m_res_arg+1),beta_arg)
    fac3=-fac*pper_arg*EwR*jn((m_res_arg-1),beta_arg)
    
    tmp=(fac1+fac2+fac3)*np.sin(eta_arg)
    return tmp

#### non linear effects####

def nonlinear_S(H,wtsq):
    tmp=H/wtsq
    return tmp

def nonlinear_H(pper,ppar,kpar,gamma,m_res,m,wce,dkpar_dt,dwcdz,dwdt):
    dwc_dt=(ppar/(gamma*m))*dwcdz
    fac1=(m_res/gamma)*dwc_dt
    fac2=-(ppar/(gamma*m))*dkpar_dt
    fac3=((kpar*pper*pper)/(2*gamma*gamma*m*m*wce))*dwcdz
    tmp=fac1+fac2+fac3-dwdt
    return tmp

def nonlinear_theta(C0,Cp1,Cm1,m_res,beta):
    tmp=C0*jn((m_res),beta)+Cp1*jn((m_res+1),beta)+Cm1*jn((m_res-1),beta)
    tmpwtsq=np.abs(tmp)
    return tmp,tmpwtsq

def nonlinear_C0(ppar,m_res,wce,kz,gamma,Ezw):
    tau=m_res*wce-(kz*ppar/const.me)
    facres=(-1)**(m_res-1)
    fac1=const.qe*kz/(gamma*const.me)
    fac2=(tau*const.qe*ppar)/(gamma*gamma*gamma*const.me*const.me*const.c_light*const.c_light)
    tmp=-facres*(fac1+fac2)*Ezw
    
    return tmp

def nonlinear_C1p(pper_arg,ppar_arg,w2_arg,Exw_arg,Eyw_arg,m_res_arg,wce_arg,kz_arg,gamma_arg):
    EwL=0.5*(Exw_arg-Eyw_arg)
    tau=m_res_arg*wce_arg-(kz_arg*ppar_arg/const.me)
    fac_mres=(-1)**(m_res_arg-1)

    fac1a=(tau*const.qe*pper_arg*EwL)/(gamma_arg*gamma_arg*gamma_arg*const.me*const.me*const.c_light*const.c_light)
    fac1b=(pper_arg*kz_arg*w2_arg)/(gamma_arg*gamma_arg*const.me)

    tmp=fac_mres*(fac1a+fac1b)
    return tmp

def nonlinear_C1m(pper_arg,ppar_arg,w1_arg,Exw_arg,Eyw_arg,m_res_arg,wce_arg,kz_arg,gamma_arg):
    EwR=0.5*(Exw_arg+Eyw_arg)
    tau=m_res_arg*wce_arg-(kz_arg*ppar_arg/const.me)
    fac_mres=(-1)**(m_res_arg-1)

    fac1a=(tau*const.qe*pper_arg*EwR)/(gamma_arg*gamma_arg*gamma_arg*const.me*const.me*const.c_light*const.c_light)
    fac1b=(pper_arg*kz_arg*w1_arg)/(gamma_arg*gamma_arg*const.me)

    tmp=fac_mres*(fac1a-fac1b)
    return tmp

