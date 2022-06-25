"""
waveproperties_mod.refr_index_warm

**Description**:
_____________________________________________________________________________________________________________________

Calculate the refractive index accounting warm plasma corrections
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

K_e:the electron warm dielectric tensor compoments

K_H:the hydrogen warm dielectric tensor compoments

K_He:the helium warm dielectric tensor compoments

K_O: the oxygen warm dielectric tensor compoments

psi_arg: wave normal angle in rad

Te: electron temperature in eV

Ti: ion temperature in eV

K110: K11 componento of the cold dielectric tensor

K220: K22 componento of the cold dielectric tensor

K330: K33 componento of the cold dielectric tensor

K120: K12 componento of the cold dielectric tensor

w_wave_arg: wave frequency

______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

ref_ind: refractive index

kappa: warm wave number

kappa_par: warm parallel wave number

kappa_per: warm perpendicular wave number
_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Kulkarni, P., Gołkowski, M., Inan, U., and Bell, T. (2015). The effect of electron and ion temperature on
the refractive index surface of 1–10 khz whistler mode waves in the inner magnetosphere. Journal of
Geophysical Research: Space Physics 120, 581–591
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""
import numpy as np

from WPIT.Environment_mod import const



def refr_index_warm(K_e,K_H,K_He,K_O,psi_arg,Te,Ti,K110,K220,K330,K120,w_wave_arg):
    
    Te_kelvin=Te*11604.5250061598  
    Ti_kelvin=Ti*11604.5250061598  
    qeT=(const.kb*Te_kelvin)/(const.me*const.c_light*const.c_light)
    qHT=(const.kb*Ti_kelvin)/(const.mH*const.c_light*const.c_light)
    qHeT=(const.kb*Ti_kelvin)/(const.mHe*const.c_light*const.c_light)
    qOT=(const.kb*Ti_kelvin)/(const.mO*const.c_light*const.c_light)
    
    K11e=K_e[0]
    K12e=K_e[1]
    K13e=K_e[2]
    K21e=K_e[3]
    K22e=K_e[4]
    K23e=K_e[5]
    K31e=K_e[6]
    K32e=K_e[7]
    K33e=K_e[8]
    
    A1e=K11e*np.sin(psi_arg)*np.sin(psi_arg)+K33e*np.cos(psi_arg)*np.cos(psi_arg)+2*K13e*np.sin(psi_arg)*np.cos(psi_arg)
    B1efac1=-(K11e*K220+K22e*K110+2*K120*K12e)*np.sin(psi_arg)*np.sin(psi_arg)
    B1efac2=-K33e*(K110+K220*np.cos(psi_arg)*np.cos(psi_arg))
    B1efac3=-K330*(K11e+K22e*np.cos(psi_arg)*np.cos(psi_arg))
    B1efac4=2*np.sin(psi_arg)*np.cos(psi_arg)*(K120*K23e-K13e*K220)
    B1e=B1efac1+B1efac2+B1efac3+B1efac4
    C1e=K33e*(K120*K120+K110*K220)+K330*(2*K120*K12e+K110*K22e+K11e*K220)
    
    K11H=K_H[0]
    K12H=K_H[1]
    K13H=K_H[2]
    K21H=K_H[3]
    K22H=K_H[4]
    K23H=K_H[5]
    K31H=K_H[6]
    K32H=K_H[7]
    K33H=K_H[8]
    
    A1H=K11H*np.sin(psi_arg)*np.sin(psi_arg)+K33H*np.cos(psi_arg)*np.cos(psi_arg)+2*K13H*np.sin(psi_arg)*np.cos(psi_arg)
    B1Hfac1=-(K11H*K220+K22H*K110+2*K120*K12H)*np.sin(psi_arg)*np.sin(psi_arg)
    B1Hfac2=-K33H*(K110+K220*np.cos(psi_arg)*np.cos(psi_arg))
    B1Hfac3=-K330*(K11H+K22H*np.cos(psi_arg)*np.cos(psi_arg))
    B1Hfac4=2*np.sin(psi_arg)*np.cos(psi_arg)*(K120*K23H-K13H*K220)
    B1H=B1Hfac1+B1Hfac2+B1Hfac3+B1Hfac4
    C1H=K33H*(K120*K120+K110*K220)+K330*(2*K120*K12H+K110*K22H+K11H*K220)
    
    
    K11He=K_He[0]
    K12He=K_He[1]
    K13He=K_He[2]
    K21He=K_He[3]
    K22He=K_He[4]
    K23He=K_He[5]
    K31He=K_He[6]
    K32He=K_He[7]
    K33He=K_He[8]
    
    A1He=K11He*np.sin(psi_arg)*np.sin(psi_arg)+K33He*np.cos(psi_arg)*np.cos(psi_arg)+2*K13He*np.sin(psi_arg)*np.cos(psi_arg)
    B1Hefac1=-(K11He*K220+K22He*K110+2*K120*K12He)*np.sin(psi_arg)*np.sin(psi_arg)
    B1Hefac2=-K33He*(K110+K220*np.cos(psi_arg)*np.cos(psi_arg))
    B1Hefac3=-K330*(K11He+K22He*np.cos(psi_arg)*np.cos(psi_arg))
    B1Hefac4=2*np.sin(psi_arg)*np.cos(psi_arg)*(K120*K23He-K13He*K220)
    B1He=B1Hefac1+B1Hefac2+B1Hefac3+B1Hefac4
    C1He=K33He*(K120*K120+K110*K220)+K330*(2*K120*K12He+K110*K22He+K11He*K220)    
    
    K11O=K_O[0]
    K12O=K_O[1]
    K13O=K_O[2]
    K21O=K_O[3]
    K22O=K_O[4]
    K23O=K_O[5]
    K31O=K_O[6]
    K32O=K_O[7]
    K33O=K_O[8]
    
    A1O=K11O*np.sin(psi_arg)*np.sin(psi_arg)+K33O*np.cos(psi_arg)*np.cos(psi_arg)+2*K13O*np.sin(psi_arg)*np.cos(psi_arg)
    B1Ofac1=-(K11O*K220+K22O*K110+2*K120*K12O)*np.sin(psi_arg)*np.sin(psi_arg)
    B1Ofac2=-K33O*(K110+K220*np.cos(psi_arg)*np.cos(psi_arg))
    B1Ofac3=-K330*(K11O+K22O*np.cos(psi_arg)*np.cos(psi_arg))
    B1Ofac4=2*np.sin(psi_arg)*np.cos(psi_arg)*(K120*K23O-K13O*K220)
    B1O=B1Ofac1+B1Ofac2+B1Ofac3+B1Ofac4
    C1O=K33O*(K120*K120+K110*K220)+K330*(2*K120*K12O+K110*K22O+K11O*K220)

    A0=K110*np.sin(psi_arg)*np.sin(psi_arg)+K330*np.cos(psi_arg)*np.cos(psi_arg)
    B0=-(K110*K220+K120*K120)*np.sin(psi_arg)*np.sin(psi_arg)-K330*(K110+K220*np.cos(psi_arg)*np.cos(psi_arg))
    C0=K330*(K120*K120+K110*K220)
    
    Awarm=(qeT*A1e+qHT*A1H+qHeT*A1He+qOT*A1O)
    Bwarm=(A0+qeT*B1e+qHT*B1H+qHeT*B1He+qOT*B1O)
    Cwarm=(B0+qeT*C1e+qHT*C1H+qHeT*C1He+qOT*C1O)
    
    # print(Awarm, Bwarm, Cwarm, C0)
        
    musq_numpy=np.roots([Awarm, Bwarm, Cwarm, C0])
    pos = np.where(musq_numpy > 0)
    ref_ind=np.sqrt(musq_numpy[pos])
    
    ref_ind=np.real(ref_ind)[0]
    
    kappa=ref_ind*w_wave_arg/const.c_light
    kappa_par=kappa*np.cos(psi_arg)
    kappa_per=kappa*np.sin(psi_arg)
    
    
    return ref_ind,kappa,kappa_par,kappa_per