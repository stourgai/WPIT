import numpy as np
from scipy.special import jn 
import matplotlib.pyplot as plt
import os
import sys
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "../../../../environment")
sys.path.append(fpath)
# print(current_dir)
# print(fpath)
import  const

def Bell_params(m_res_arg,ppar_arg,pper_arg,Bxw_arg, Byw_arg,Exw_arg,Eyw_arg,Ezw_arg,kz_arg,kx_arg,wce_arg):
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


def z_rk_bort(ppar_arg,gamma_arg):
    krk=ppar_arg/(gamma_arg*const.me)
    return krk

def p_par_rk_bort(pper_arg,eta_arg,wtau_sq_arg,kz_arg,gamma_arg,wce_arg,dwds_arg):
    lrk=((wtau_sq_arg*const.me)/kz_arg)*np.sin(eta_arg)-(1/(const.me*gamma_arg))*((pper_arg*pper_arg)/(2*wce_arg))*dwds_arg
    return lrk    

def p_per_rk_bort(ppar_arg,pper_arg,eta_arg,w1_arg,w2_arg,beta_arg,gamma_arg,R1_arg,R2_arg,m_res_arg,wce_arg,dwds_arg):
    mrk=(-(pow((-1),(m_res_arg-1))*(w1_arg*((ppar_arg/gamma_arg)+(const.me*R1_arg))*jn((m_res_arg-1),beta_arg)-
            w2_arg*((ppar_arg/gamma_arg)-(const.me*R2_arg))*jn((m_res_arg+1),beta_arg))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg*ppar_arg)/(2*wce_arg))*dwds_arg)
    return mrk

def eta_rk_bort(ppar_arg,m_res_arg,wce_arg,wwave_arg,gamma_arg,kz_arg):
    nrk=((m_res_arg*wce_arg)/gamma_arg)-wwave_arg-kz_arg*(ppar_arg/(const.me*gamma_arg))
    return nrk

def lamda_rk_bort(ppar_arg,lamda_arg,gamma_arg,L_arg):
    ork=ppar_arg/(gamma_arg*const.me*L_arg*const.Re*np.sqrt(1+3*np.sin(lamda_arg)*
                    np.sin(lamda_arg))*np.cos(lamda_arg))
    return ork

def alpha_rk_bort(pper_arg,eta_arg,alpha_arg,wtau_sq_arg,kz_arg,m_res_arg,wce_arg,w_wave_arg,gamma_arg,dwhds_arg):

    prk=(-((const.me*wtau_sq_arg/(kz_arg*pper_arg))*(1+((np.cos(alpha_arg)*np.cos(alpha_arg))/(m_res_arg*(wce_arg/w_wave_arg)-1))))*np.sin(eta_arg)+
            ((1/(const.me*gamma_arg))*(pper_arg/(2*wce_arg)))*dwhds_arg)

    return prk

def alpha_eq_rk_bort(ppar_arg,pper_arg,alpha_arg,aeq_arg,eta_arg,w1_arg,R1_arg,w2_arg,R2_arg,gamma_arg,beta_arg,wtausq_arg,kz_arg,m_res_arg):
    pmag=np.sqrt(ppar_arg**2+pper_arg**2)
    fac1=-(1/(pmag*pmag))*(np.tan(aeq_arg)/np.tan(alpha_arg))
    fac2a=w1_arg*((ppar_arg/gamma_arg)+const.me*R1_arg)*jn((m_res_arg-1),beta_arg)
    fac2b=w2_arg*((ppar_arg/gamma_arg)-const.me*R2_arg)*jn((m_res_arg+1),beta_arg)
    fac2=(fac2a-fac2b)*ppar_arg
    fac3=wtausq_arg*const.me*pper_arg/kz_arg
    tmp=fac1*(fac2+fac3)*np.sin(eta_arg)
    
    return tmp

