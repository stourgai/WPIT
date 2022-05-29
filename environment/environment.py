import const 
import numpy as np

def B_dipole(L_arg,lamda_arg,r_arg):
    
    B0=31200*10**(-9)
    Br_tmp=-2*B0*((const.Re/r_arg)**3)*np.sin(lamda_arg)
    Bl_tmp=B0*((const.Re/r_arg)**3)*np.cos(lamda_arg)
    Bt_tmp=0
    Bmag_tmp=np.sqrt(Br_tmp*Br_tmp+Bl_tmp*Bl_tmp+Bt_tmp*Bt_tmp)
    return Br_tmp,Bl_tmp,Bt_tmp,Bmag_tmp

def Bmag_dipole(L_arg,lamda_arg):
    #----calculate the dipole magnetic field strength
    #L_arg geomagnetic L shell
    #lamda_arg geomagnetic latitude in rad
    
    B0=31200*10**(-9)
    slat=np.sin(lamda_arg)
    clat=np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    Bmag=(B0/(L_arg**3))*slat_term/(clat**6)

    return Bmag

def cyclotron(B_arg,q_arg,m_arg):
    #----calculate the cyclotron frequency
    #q_arg the species charge in Coulomb
    #m_arg the species mass in kg
    #B_arg the geomagnetic field strength in Tesla
    omega_tmp=(np.abs(q_arg)*B_arg)/m_arg
    fc_tmp=omega_tmp/(2*np.pi)
    return omega_tmp, fc_tmp

def dB_ds(B_arg,lamda_arg,L_arg):

    slat = np.sin(lamda_arg)
    clat = np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    dB_ds_arg = (3.0*B_arg/(L_arg*const.Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))

    return dB_ds_arg

def densities_denton(ne0_arg,lamda_arg):
    #-----calculate species densities (assuming 94%H+, 5.4%He+, 0.6%O+)
    #lamda geomagnetic latitude
    clat=np.cos(lamda_arg)
    n_e_tmp=ne0_arg*(clat**(-4))

    return n_e_tmp

def dipole_cart(xtmp,ytmp,ztmp):
    Re=6378137 #Earth mean radius in meters
    mu=np.pi*4*10**-7 #permiability of free space in H/m
    B0=3.12*10**-5 #reference magnetic field in T
    M=(4*np.pi*B0*Re**3)/mu #dipole moment in Am^2

    rtmp=np.sqrt(xtmp**2+ytmp**2+ztmp**2)
    theta_tmp=np.arccos(ztmp/rtmp)

    if xtmp>0:
        phi_tmp=np.arctan(ytmp/xtmp)
    elif xtmp<0 and ytmp>=0:
        phi_tmp=np.arctan(ytmp/xtmp)+np.pi 
    elif xtmp<0 and ytmp<0:
        phi_tmp=np.arctan(ytmp/xtmp)-np.pi 
    elif xtmp == 0 and ytmp>0:
        phi_tmp=np.pi/2
    elif xtmp == 0 and ytmp<0:
        phi_tmp=-np.pi/2
    elif xtmp == 0  and ytmp == 0:
        phi_tmp =0

    rtmp_unit=[np.sin(theta_tmp)*np.cos(phi_tmp),np.sin(theta_tmp)*np.sin(phi_tmp),np.cos(theta_tmp)]
    theta_tmp_unit=[np.cos(theta_tmp)*np.cos(phi_tmp),np.cos(theta_tmp)*np.sin(phi_tmp),-np.sin(theta_tmp)]
    phi_tmp_unit=[-np.sin(phi_tmp),np.cos(phi_tmp),0.0]

    Br_mag=((-2*mu*M*np.cos(theta_tmp)/(4*np.pi*rtmp**3)))
    Btheta_mag=((-mu*M*np.sin(theta_tmp)/(4*np.pi*rtmp**3)))
    Bphi_mag=0.0

    Br=Br_mag*np.asarray(rtmp_unit)
    Btheta=Btheta_mag*np.asarray(theta_tmp_unit)
    Bphi=Bphi_mag*np.asarray(phi_tmp_unit)

    BX=Br[0]+Btheta[0]+Bphi[0] 
    BY=Br[1]+Btheta[1]+Bphi[1] 
    BZ=Br[2]+Btheta[2]+Bphi[2]

    Bmag=np.sqrt(BX**2+BY**2+BZ**2) 
    bunit=[BX/Bmag,BY/Bmag,BZ/Bmag]
    # print (BX,BY,BZ,Bmag)
    return BX, BY, BZ, Bmag, Br_mag, Btheta_mag 

def dwc_ds(wc_arg,lamda_arg,L_arg):

    slat = np.sin(lamda_arg)
    clat = np.cos(lamda_arg)
    slat_term = np.sqrt(1 + 3*slat*slat)
    dwce_ds_arg = (3.0*wc_arg/(L_arg*const.Re)) *(slat/slat_term) * (1.0/(slat_term*slat_term) + 2.0/(clat*clat))
    return dwce_ds_arg

def f_lower_hybrid(fc_arg,fp_arg,ne_arg,nH_arg,nO_arg,nHe_arg):
    ion_fac=(((nH_arg)/const.mH)+((nHe_arg)/const.mHe)+((nO_arg)/const.mO))
    m_fac=(const.me/(ne_arg))
    tmp=np.sqrt(m_fac*((fc_arg*fc_arg*fp_arg*fp_arg)/(fc_arg*fc_arg+fp_arg*fp_arg)))
    return tmp

def geo_lat2geod_lat(phi):
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geod_lat = np.arctan(np.tan(phi_rad) / (1 - e ** 2))
    geod_lat = np.rad2deg(geod_lat)
    return geod_lat  # in degrees

def geo2geod(lat_geo_phi, lon_geo_lmd, alt_geo):
    lat_geod_phi = geo_lat2geod_lat(lat_geo_phi)  # degrees
    lon_geod_lmd = lon_geo_lmd  # degrees
    alt_geod = alt_geo  # km
    return lat_geod_phi, lon_geod_lmd, alt_geod

def larmor(uperp_arg,gamma_arg,B_arg):
    rl_tmp=gamma_arg*const.me*uperp_arg/(const.qe*B_arg)
    return rl_tmp

def magLshell(r_arg,lat_arg):
#     r=np.sqrt(posx[i]*posx[i]+posy[i]*posy[i]+posz[i]*posz[i])
    L_tmp=r_arg/(const.Re*np.cos(lat_arg)*np.cos(lat_arg))
    return L_tmp

def omega_lower_hybrid(wce_arg,wpe_arg,wci_arg,wpi_arg):
    fac1=(wce_arg*wce_arg+wpe_arg*wpe_arg+wci_arg*wci_arg+wpi_arg*wpi_arg)/2
    fac2=(wce_arg*wce_arg+wpe_arg*wpe_arg-wci_arg*wci_arg-wpi_arg*wpi_arg)**2+4*wpe_arg*wpe_arg*wpi_arg*wpi_arg
    tmp=fac1-0.5*fac2
    return tmp

def omega_lower_hybrid_v2(wce_arg,wpe_arg,wci_arg,wpi_arg):
    fac1=1/(wci_arg*wce_arg)
    fac2=1/(wpi_arg**2)
    sqrtfac=np.sqrt(fac1+fac2)
    tmp=1/sqrtfac
    return tmp

def omega_plasma(n_arg,q_arg,m_arg):
    omegap_tmp=np.sqrt((n_arg*q_arg*q_arg)/(m_arg*const.epsilon0))
    fp_tmp=omegap_tmp/(2*np.pi)
    return omegap_tmp,fp_tmp

def omega_upper_hybrid(wce_arg,wpe_arg,wci_arg,wpi_arg):
    fac1=(wce_arg*wce_arg+wpe_arg*wpe_arg+wci_arg*wci_arg+wpi_arg*wpi_arg)/2
    fac2=(wce_arg*wce_arg+wpe_arg*wpe_arg-wci_arg*wci_arg-wpi_arg*wpi_arg)**2+4*wpe_arg*wpe_arg*wpi_arg*wpi_arg
    wuh_tmp=fac1+0.5*fac2
    return wuh_tmp

