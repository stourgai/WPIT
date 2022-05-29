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
