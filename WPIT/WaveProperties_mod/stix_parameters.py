"""
waveproperties_mod.stix_parameters

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the the Stix parameters

_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

w_arg: wave frequency in Hz

Ne_arg: electron number density in m^-3

NH_arg: H+ number density in m^-3 (if available, else 0)

NHe_arg: He+ number density in m^-3 (if available, else 0)

NO_arg: O+ number density in m^-3 (if available, else 0)

B0mag_arg: Geomagnetic field in T
_______________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

Sstix: S parameter

Dstix: D parameter

Pstix: P parameter

Rstix: R parameter

Lstix: L parameter
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Stix, T. H. (1992). Waves in plasmas (Springer Science & Business Media)
_________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

from WPIT.Environment_mod import const


def stix_parameters(w_arg, Ne_arg, NH_arg, NHe_arg, NO_arg, B0mag_arg):

    wpse_arg = (Ne_arg*(const.qe*const.qe)/(const.me*const.epsilon0))  #electron plasma frequency squared
    wce_arg = (const.qe*B0mag_arg)/const.me        #electron cyclotron frequency

    wpsH_arg = (NH_arg*(const.qi*const.qi)/(const.mH*const.epsilon0))  #hydrogen plasma frequency squared
    wcH_arg = (const.qi*B0mag_arg)/const.mH       #hydrogen cyclotron frequency 

    wpsHe_arg = (NHe_arg*(const.qi*const.qi)/(const.mHe*const.epsilon0))   #helium plasma frequency squared
    wcHe_arg = (const.qi*B0mag_arg)/const.mHe         #helium cyclotron frequency 

    wpsO_arg = (NO_arg*(const.qi*const.qi)/(const.mO*const.eps0))   #oxygen plasma frequency squared
    wcO_arg = (const.qi*B0mag_arg)/const.mO        #oxugen cyclotron frequency 
    
    Rface=(wpse_arg/(w_arg*(w_arg-wce_arg)))
    RfacH=(wpsH_arg/(w_arg*(w_arg+wcH_arg)))
    RfacHe=(wpsHe_arg/(w_arg*(w_arg+wcHe_arg)))
    RfacO=(wpsO_arg/(w_arg*(w_arg+wcO_arg)))

    Lface=(wpse_arg/(w_arg*(w_arg+wce_arg)))
    LfacH=(wpsH_arg/(w_arg*(w_arg-wcH_arg)))
    LfacHe=(wpsHe_arg/(w_arg*(w_arg-wcHe_arg)))
    LfacO=(wpsO_arg/(w_arg*(w_arg-wcO_arg)))
    
    Rstix=1-Rface-RfacH-RfacO-RfacHe
    Lstix=1-Lface-LfacH-LfacO-LfacHe
    Pstix=1-(wpse_arg/w_arg**2)-(wpsH_arg/w_arg**2)-(wpsHe_arg/w_arg**2)-(wpsO_arg/w_arg**2)
    Sstix = (Rstix+Lstix)/2
    Dstix = (Rstix-Lstix)/2

    return Sstix,Dstix,Pstix,Rstix,Lstix