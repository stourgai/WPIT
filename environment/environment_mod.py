#########


########


import numpy as np
import const

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
  
 
def geo2geod(lat_geo_phi, lon_geo_lmd, alt_geo):
    lat_geod_phi = geo_lat2geod_lat(lat_geo_phi)  # degrees
    lon_geod_lmd = lon_geo_lmd  # degrees
    alt_geod = alt_geo  # km
    return lat_geod_phi, lon_geod_lmd, alt_geod
  
def geo_lat2geod_lat(phi):
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geod_lat = np.arctan(np.tan(phi_rad) / (1 - e ** 2))
    geod_lat = np.rad2deg(geod_lat)
    return geod_lat  # in degrees
  
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

def aeq2alpha(L_arg,lambda_arg,aeq_arg):
    Blam0=Bmag_dipole(L_arg,lambda_arg)
#    print(Blam0)
    Beq0=Bmag_dipole(L_arg,0)
#    print(Beq0)
    salpha0=np.sin(aeq_arg)*np.sqrt(Blam0/Beq0)
#     print(np.rad2deg(salpha0))
    alpha0=np.arcsin(salpha0)
    
    return alpha0

def alpha2aeq(L_arg,lambda_arg,alpha_arg):
    Blam0=Bmag_dipole(L_arg,lambda_arg)
    Beq0=Bmag_dipole(L_arg,0)
    salphaeq0=np.sin(alpha_arg)*np.sqrt(Beq0/Blam0)

    alphaeq0=np.arcsin(salphaeq0)
    return alphaeq0


def momentums(Ekev,alpha,m_arg):
    Ejoule0=1.602176487E-16*Ekev
    gamma0=(Ejoule0/(m_arg*(const.c_light**2))) +1
    speed0=np.sqrt(1- (1/(gamma0**2)))*const.c_light
    upar0=speed0*np.cos(alpha)
    uper0=speed0*np.sin(alpha)
    pper0=gamma0*m_arg*uper0
    ppar0=gamma0*m_arg*upar0
    
    return upar0,uper0,ppar0,pper0,gamma0

import numpy as np
import matplotlib.pyplot as plt
def L_pp(Kpmax):
    tmp=5.6-0.46*Kpmax
    return tmp

#saturated plasmasphere 2.25<=L<=Lppi
def sat_plasmasphere(L,d,R):
    log_ne=(-0.3145*L+3.9043)+(0.15*(np.cos(2*np.pi*(d+9)/365)-0.5*np.cos(4*np.pi*(d+9)/365))+0.00127*R-0.0635)*np.exp(-(L-2/1.5))
    ne=10**log_ne
#     print(ne)
    return ne

#Lppi<=L<=Lppo
def plasmapause(L,Lppi,ne_lppi,mlt,d,R):
    if mlt>=0 and mlt<6:
        ne=ne_lppi*10**(-(L-Lppi)/0.1)
    else:
        ne=ne_lppi*10**(-(L-Lppi)/(0.1+0.011*(mlt-6)))
    return ne

#extended plasma trough 2.25<=L<=8
def ext_trough(L,mlt):
    if mlt>=0 and mlt<6:
        ne=(5800+300*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    else:
        ne=(-800+1400*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    return ne  

def trough(ne_Lppo,L,Lppo):
    ne=ne_Lppo*((L/Lppo)**(-4.5))+(1-np.exp(-(L-2)/10))
    return ne

#find the Lppo from the cross point of the two functions
def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc





def carpender_anderson(Lsh,Kpmax,day,mlt,Rb):

    L_plasma=[]
    ne_plasma=[]
    ne_trough=[]
    ne_final=[]
    L_final=[]
    # ne_final.append(10**6)
    # L_final.append(1)
    L_array=np.arange(2.25,8,0.001)

    L_ppi=L_pp(Kpmax)
    ne_lppi=sat_plasmasphere(L_ppi,day,Rb)

    for i in range(0,len(L_array)):
        if L_array[i]<L_ppi:
            ne=sat_plasmasphere(L_array[i],day,Rb)
            ne_plasma.append(ne)
            L_plasma.append(L_array[i])
        if L_array[i]>L_ppi:
            ne=plasmapause(L_array[i],L_ppi,ne_lppi,mlt,day,Rb)
            ne_plasma.append(ne)
            L_plasma.append(L_array[i])
    for i in range(0,len(L_array)):
        ne_tr=ext_trough(L_array[i],mlt)
        ne_trough.append(ne_tr)

    ne_plasma=np.asarray(ne_plasma)
    ne_trough=np.asarray(ne_trough)

    xc, yc = interpolated_intercept(L_array,ne_plasma,ne_trough)
    L_ppo=xc[0][0]

    ne_Lppo=plasmapause(L_ppo,L_ppi,ne_lppi,mlt,day,Rb)
    # print(ne_Lppo)
#     for i in range(0,len(L_array)):
#         if L_array[i]<L_ppi:
#             ne=sat_plasmasphere(L_array[i],day,Rb)
#             ne_final.append(ne)
#             L_final.append(L_array[i])
#         if L_ppi<=L_array[i]<=L_ppo:
#             ne=plasmapause(L_array[i],L_ppi,ne_lppi,mlt,day,Rb)
#             ne_final.append(ne)
#             L_final.append(L_array[i])
#         if L_ppo<=L_array[i]<=8:
#             ne=trough(ne_Lppo,L_array[i],L_ppo)
#             ne_final.append(ne)
#             L_final.append(L_array[i]) 
    
    if Lsh<L_ppi:
        ne_eq=sat_plasmasphere(Lsh,day,Rb)
    if L_ppi<=Lsh<=L_ppo:
        ne_eq=plasmapause(Lsh,L_ppi,ne_lppi,mlt,day,Rb)
    if L_ppo<=Lsh<=8:
        ne_eq=trough(ne_Lppo,Lsh,L_ppo)
        
    return ne_eq

def densities_ozhogin(L_arg,lambda_arg):
    neqtmp=10**(4.4693-0.4903*L_arg)
    lamda_inv=np.arccos(1/L_arg)
    n_lamtmp=neqtmp*((np.cos((np.pi/2)*(lambda_arg/lamda_inv)))**(-0.75))
    return neqtmp,n_lamtmp
