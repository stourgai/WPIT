"""
initial_velocity.dB_ds

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate velocities and momentums from energy and pitch angle

_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

Ekev: energy in keV

alpha: local pitch angle in rad

L_arg: L shell
_______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

upar0: parallel velocity in m/s

uper0: perpendicular velocity in m/s

ppar0: parallel momentum in N*s

pper0: perpendicular momentum in N*s

gamma0: Lorentz factor
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np
from environment_mod import const

def initial_velocity(Ekev,alpha,m_arg):
    Ejoule0=1.602176487E-16*Ekev
    gamma0=(Ejoule0/(m_arg*(const.c_light**2))) +1
    speed0=np.sqrt(1- (1/(gamma0**2)))*const.c_light
    upar0=speed0*np.cos(alpha)
    uper0=speed0*np.sin(alpha)
    pper0=gamma0*m_arg*uper0
    ppar0=gamma0*m_arg*upar0
    
    return upar0,uper0,ppar0,pper0,gamma0