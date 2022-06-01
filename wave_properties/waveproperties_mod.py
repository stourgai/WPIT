import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/../environment")
sys.path.append(fpath)

import  const

def stix_parameters(w_arg, Ne_arg, NH_arg, NHe_arg, NO_arg, B0mag_arg):
    #----calculate Stix paramet_arg
    #N*_arg  densities in m^-3
    #w_arg wave frequency in Hz
    #B0mag magnitude of the ambient magnetic field in Tesla

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
    
    R_arg=1-Rface-RfacH-RfacO-RfacHe
    L_arg=1-Lface-LfacH-LfacO-LfacHe
    P_arg=1-(wpse_arg/w_arg**2)-(wpsH_arg/w_arg**2)-(wpsHe_arg/w_arg**2)-(wpsO_arg/w_arg**2)
    S_arg = (R_arg+L_arg)/2
    D_arg = (R_arg-L_arg)/2

    return S_arg,D_arg,P_arg,R_arg,L_arg
    
    
def ref_index(theta,S,P,R,L):
    A=S*np.sin(theta)*np.sin(theta)+P*np.cos(theta)*np.cos(theta)
    B=R*L*np.sin(theta)*np.sin(theta)+P*S*(1+np.cos(theta)*np.cos(theta))
    C=P*R*L
    F=np.sqrt(B*B-4*A*C)
    eta_sq_plus=(B+F)/(2*A)
    eta_sq_minus=(B-F)/(2*A)
    return eta_sq_plus,eta_sq_minus

def dispersion_appleton(w_arg,wpe_arg,wce_arg,theta_arg):
    fac1=(wpe_arg*wpe_arg)/(w_arg*w_arg)
    fac2=(wce_arg*wce_arg*np.sin(theta_arg)*np.sin(theta_arg))/(2*(w_arg*w_arg-wpe_arg*wpe_arg))
    sqrtfac=np.sqrt(fac2*fac2+(((wce_arg*wce_arg)/(w_arg*w_arg))*np.cos(theta_arg)*np.cos(theta_arg)))
    denomfac_plus=1-fac2+sqrtfac
    denomfac_minus=1-fac2-sqrtfac
    musq_plus=1-(fac1/denomfac_plus)
    musq_minus=1-(fac1/denomfac_minus)
    
    return musq_plus,musq_minus  
    
    
def whislter_amplitudes_bell(mu_arg,P_arg,D_arg,S_arg,Byw_arg,theta_arg):

    mu_sq_arg=mu_arg**2
    fac1= (P_arg-mu_sq_arg*(np.sin(theta_arg)**2)) 
    Byw_arg=Byw_arg

    Bxw_arg=(-(D_arg*fac1)/(P_arg*(S_arg-mu_arg**2)))*Byw_arg
    Bzw_arg=((D_arg*np.sin(theta_arg)*fac1)/(P_arg*np.cos(theta_arg)*(S_arg-mu_arg**2)))*Byw_arg
    Exw_arg=((const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg))*Byw_arg)
    Eyw_arg=((D_arg*const.c_light*fac1)/(mu_arg*P_arg*np.cos(theta_arg)*(mu_arg**2-S_arg)))*Byw_arg
    Ezw_arg=(-(const.c_light*mu_arg*np.sin(theta_arg))/P_arg)*Byw_arg
    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg

def whistler_waves_bell(Bxwc_arg, Bywc_arg, Bzwc_arg, Exwc_arg, Eywc_arg, Ezwc_arg,time_arg,kz_arg,zeta_arg,kx_arg,chi_arg,w_wave_arg):
    Phi_arg=(w_wave_arg*time_arg-kz_arg*zeta_arg-kx_arg*chi_arg)   #should i add kx*x ?

    Exw_arg=-Exwc_arg*np.sin(Phi_arg)  
    Eyw_arg=Eywc_arg*np.cos(Phi_arg)
    Ezw_arg=-Ezwc_arg*np.sin(Phi_arg)
    Bxw_arg=Bxwc_arg*np.cos(Phi_arg)  
    Byw_arg=Bywc_arg*np.sin(Phi_arg)
    Bzw_arg=-Bzwc_arg*np.cos(Phi_arg)

    return Bxw_arg, Byw_arg, Bzw_arg, Exw_arg, Eyw_arg, Ezw_arg, Phi_arg

def wave_packet(Bw0_arg,lamda_arg):
    Bwave=Bw0_arg*(np.tanh(-2*np.rad2deg(lamda_arg)-1)+1)/2
#     print((np.tanh(-2*np.rad2deg(lamda_arg)-1)+1)/2)
    return Bwave

def whislter_amplitudes_li(mu,P,D,S,Bw_tot_li,psi):
    #Li uses a different approach in defining the wave fields (see notebook for more details)

    I_w=(Bw_tot_li/(mu*np.sqrt(D*D*(P-mu*mu*np.sin(psi)*np.sin(psi))*(P-mu*mu*np.sin(psi)*np.sin(psi))
                           +P*P*np.cos(psi)*np.cos(psi)*(S-mu*mu)*(S-mu*mu))))
    mu_sq_li=mu*mu
    fac1= (P-mu_sq_li*np.sin(psi)*np.sin(psi)) 
    fac2=S-mu*mu
    Exw_li=const.c_light*I_w*fac1*fac2
    Eyw_li=const.c_light*I_w*D*fac1
    Ezw_li=-const.c_light*I_w*mu_sq_li*np.cos(psi)*np.sin(psi)*fac2
    Bxw_li=-I_w*D*np.cos(psi)*fac1*mu
    Byw_li=I_w*P*np.cos(psi)*fac2*mu
    Bzw_li=I_w*D*np.sin(psi)*fac1*mu
                
    return Bxw_li, Byw_li, Bzw_li, Exw_li, Eyw_li, Ezw_li

def whistler_waves_li(Bxwc, Bywc, Bzwc, Exwc, Eywc, Ezwc,time,kz,kx,zeta,chi,w_wave):
    Phi=(w_wave*time-kz*zeta-kx*chi) 

    Exw=Exwc*np.sin(Phi)  
    Eyw=Eywc*np.cos(Phi)
    Ezw=Ezwc*np.sin(Phi)
    Bxw=Bxwc*np.cos(Phi)  
    Byw=Bywc*np.sin(Phi)
    Bzw=Bzwc*np.cos(Phi)

    return Bxw, Byw, Bzw, Exw, Eyw, Ezw, Phi

def res_angle(P_arg,S_arg):
    tansq=-P_arg/S_arg
    tan=np.sqrt(tansq)
    thetares=np.arctan(tan)
    return thetares

def resonant_velocity(m_res_arg,w_wave_arg,kz_arg,wce_arg,alpha_arg):
    fac1 = w_wave_arg*w_wave_arg*kz_arg*kz_arg
    fac2 = ((m_res_arg*wce_arg)*(m_res_arg*wce_arg))-(w_wave_arg*w_wave_arg)
    fac3 = kz_arg*kz_arg + ((m_res_arg*wce_arg)*(m_res_arg*wce_arg)) /((const.c_light*np.cos(alpha_arg))*(const.c_light*np.cos(alpha_arg)))

    if m_res_arg == 0:
        direction=-1.0*np.sign(kz_arg)
    else:
        direction = np.sign(kz_arg)*np.sign(m_res_arg)

    v_para_res = ( direction*np.sqrt(fac1 + fac2*fac3) - w_wave_arg*kz_arg) / fac3
    
    v_tot_res = v_para_res / np.cos(alpha_arg)
    v_per_res=v_tot_res*np.sin(alpha_arg)
    gamma_res=1.0/np.sqrt( 1-((v_tot_res*v_tot_res)/(const.c_light*const.c_light)))
    E_rest=const.me*const.c_light*const.c_light

    E_res=(gamma_res-1)*E_rest
    return v_para_res, v_per_res, v_tot_res, E_res,gamma_res

