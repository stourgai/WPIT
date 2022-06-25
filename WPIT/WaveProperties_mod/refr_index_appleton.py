"""
waveproperties_mod.refr_index_appleton

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the the refractive index and the wave numbers
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w_arg: wave frequency

wpe_arg: electron plasma frequency

wce_arg: electron cyclotron frequency

theta_arg: wave normal angle in rad
______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

eta_sq_plus: the plus (+) root of the dispersion relation

eta_sq_minus: the minus (-) root of the dispersion relation

ref_ind: refractive index (for the 𝑛2>0 root)

kappa: wave number

kappa_par: parallel component of the wave number

kappa_per: perpendicular component of the wave number

______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Appleton, E. V. (1932). Wireless studies of the ionosphere. Institution of Electrical Engineers-Proceedings
of the Wireless Section of the Institution 7, 257–265
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np

from WPIT.Environment_mod import const


def refr_index_appleton(w_arg,wpe_arg,wce_arg,theta_arg):
    fac1=(wpe_arg*wpe_arg)/(w_arg*w_arg)
    fac2=(wce_arg*wce_arg*np.sin(theta_arg)*np.sin(theta_arg))/(2*(w_arg*w_arg-wpe_arg*wpe_arg))
    sqrtfac=np.sqrt(fac2*fac2+(((wce_arg*wce_arg)/(w_arg*w_arg))*np.cos(theta_arg)*np.cos(theta_arg)))
    denomfac_plus=1-fac2+sqrtfac
    denomfac_minus=1-fac2-sqrtfac
    musq_plus=1-(fac1/denomfac_plus)
    musq_minus=1-(fac1/denomfac_minus)

    
    if musq_plus>0:
        ref_ind=np.sqrt(musq_plus)
    else:
        ref_ind=np.sqrt(musq_minus)

    
    kappa=ref_ind*w_arg/const.c_light
    kappa_par=kappa*np.cos(theta_arg)
    kappa_per=kappa*np.sin(theta_arg)
    return musq_plus,musq_minus,ref_ind,kappa,kappa_par,kappa_per