import numpy as np
from scipy.special import jn 

from WPIT.Environment_mod import const

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