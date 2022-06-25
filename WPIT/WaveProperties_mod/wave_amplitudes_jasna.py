"""
waveproperties_mod.wave_amplitudes_jasna

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the wave amplitudes
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

P: Stix P parameter

D: Stix D parameter

S: Stix S parameter

theta_arg: wave normal angle in rad

ref_arg: refractive index

power_arg: wave Poynting flux in mW/m^2
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

Jasna, D., Inan, U. S., and Bell, T. F. (1992). Precipitation of suprathermal (100 ev) electrons by oblique
whistler waves. Geophysical research letters 19, 1639–1642
_________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np


from WPIT.Environment_mod import const




def wave_amplitudes_jasna(P_arg,S_arg,D_arg,theta_arg,ref_arg,power_arg):
    musq=ref_arg*ref_arg
    Xstix = P_arg/(P_arg - musq*np.sin(theta_arg)*np.sin(theta_arg))
    
    rho1=((musq-S_arg)*musq*np.sin(theta_arg)*np.cos(theta_arg))/(D_arg*(musq*np.sin(theta_arg)*np.sin(theta_arg)-P_arg))
    rho2 = (musq - S_arg) / D_arg

    #print num
    Byw_sq = ( (2.0*const.mu0/const.c_light) * (power_arg*Xstix*Xstix*rho2*rho2*ref_arg*np.cos(theta_arg)) /
                np.sqrt( pow((np.tan(theta_arg)-rho1*rho2*Xstix),2) + pow((1+rho2*rho2*Xstix),2)) )


    mu_sq_arg=ref_arg**2
    fac1= (P_arg-mu_sq_arg*(np.sin(theta_arg)**2)) 
    Byw_arg=np.sqrt(Byw_sq)

    Bxw_arg=(-(D_arg*fac1)/(P_arg*(S_arg-ref_arg**2)))*Byw_arg
    Bzw_arg=((D_arg*np.sin(theta_arg)*fac1)/(P_arg*np.cos(theta_arg)*(S_arg-ref_arg**2)))*Byw_arg
    Exw_arg=((const.c_light*fac1)/(ref_arg*P_arg*np.cos(theta_arg))*Byw_arg)
    Eyw_arg=((D_arg*const.c_light*fac1)/(ref_arg*P_arg*np.cos(theta_arg)*(ref_arg**2-S_arg)))*Byw_arg
    Ezw_arg=(-(const.c_light*ref_arg*np.sin(theta_arg))/P_arg)*Byw_arg
    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg