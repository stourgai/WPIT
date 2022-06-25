"""
waveproperties_mod.wave_amplitudes_bell

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the wave amplitudes
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

mu_arg: the refractive index

P_arg: Stix P parameter

D_arg: Stix D parameter

S_arg: Stix S parameter

Byw_arg: By component of the wave in T

theta_arg: wave normal angle in rad
_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

Bxw_arg:Bx component of the wave in T

Byw_arg:By component of the wave in T

Bzw_arg:Bz component of the wave in T

Exw_arg:Ex component of the wave in V/m

Eyw_arg:Ey component of the wave in V/m

Ezw_arg:Ez component of the wave in V/m
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Bell, T. F. (1984). The nonlinear gyroresonance interaction between energetic electrons and coherent vlf
waves propagating at an arbitrary angle with respect to the earth’s magnetic field. Journal of Geophysical
Research: Space Physics 89, 905–918
_________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np

from WPIT.Environment_mod import const


def wave_amplitudes_bell(mu_arg,P_arg,D_arg,S_arg,Byw_arg,theta_arg):

    mu_sq_arg=mu_arg**2
    fac1= (P_arg-mu_sq_arg*(np.sin(theta_arg)**2)) 
    Byw_arg=Byw_arg

    Bxw_arg=(-(D_arg*fac1)/(P_arg*(S_arg-mu_arg**2)))*Byw_arg
    Bzw_arg=((D_arg*np.sin(theta_arg)*fac1)/(P_arg*np.cos(theta_arg)*(S_arg-mu_arg**2)))*Byw_arg
    Exw_arg=((const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg))*Byw_arg)
    Eyw_arg=((D_arg*const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg)*(mu_arg**2-S_arg)))*Byw_arg
    Ezw_arg=(-(const.c_light*mu_arg*np.sin(theta_arg))/P_arg)*Byw_arg
    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg