import numpy as np
from environment_mod import const
def momentums(Ekev,alpha,m_arg):
    Ejoule0=1.602176487E-16*Ekev
    gamma0=(Ejoule0/(m_arg*(const.c_light**2))) +1
    speed0=np.sqrt(1- (1/(gamma0**2)))*const.c_light
    upar0=speed0*np.cos(alpha)
    uper0=speed0*np.sin(alpha)
    pper0=gamma0*m_arg*uper0
    ppar0=gamma0*m_arg*upar0
    
    return upar0,uper0,ppar0,pper0,gamma0