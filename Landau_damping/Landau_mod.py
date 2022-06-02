import pandas as pd    
import numpy as np
from waveproperties_mod import res_angle,ref_index
import os
import sys
import datetime as dt
from scipy import integrate as scint
from spacepy import coordinates as coord
from spacepy.time import Ticktock
from scipy import special as scp
import matplotlib.pyplot as plt
current_dir =  os.path.abspath(os.path.dirname('__file__'))

fpath = os.path.abspath(current_dir + "/../environment")
sys.path.append(fpath)
# print(current_dir)
# print(fpath)
import  const

################################################
#Bell 2002

#Bell, T. F., Inan, U. S., Bortnik, J., and Scudder, J. D., 
# The Landau damping of magnetospherically reflected whistlers within the plasmasphere, 
# Geophys. Res. Lett., 29( 15), doi:10.1029/2002GL014752, 2002. 
def stix_parameters(w, qs, Ns, ms, nus, B0mag):

    # wps2 = (Ns*qs**2/ms/const.EPS0)*(w/(w))
    # wcs = ((qs*B0mag)/ms)*(w/(w))

    wps21 = (Ns[0]*qs[0]**2/ms[0]/const.eps0)*(w/(w))
    wcs1 = ((qs[0]*B0mag)/ms[0])*(w/(w)) 

    wps22 = (Ns[1]*qs[1]**2/ms[1]/const.eps0)*(w/(w))
    wcs2 = ((qs[1]*B0mag)/ms[1])*(w/(w)) 

    wps23 = (Ns[2]*qs[2]**2/ms[2]/const.eps0)*(w/(w))
    wcs3 = ((qs[2]*B0mag)/ms[2])*(w/(w)) 

    wps24 = (Ns[3]*qs[3]**2/ms[3]/const.eps0)*(w/(w))
    wcs4 = ((qs[3]*B0mag)/ms[3])*(w/(w)) 

    R=1-(wps21/(w*(w+wcs1)))-(wps22/(w*(w+wcs2)))-(wps23/(w*(w+wcs3)))-(wps24/(w*(w+wcs4)))
    L=1-(wps21/(w*(w-wcs1)))-(wps22/(w*(w-wcs2)))-(wps23/(w*(w-wcs3)))-(wps24/(w*(w-wcs4)))
    P=1-(wps21/w**2)-(wps22/w**2)-(wps23/w**2)-(wps24/w**2)
    # R = 1-sum(wps2/(w*(w+wcs)))
    # L = 1-sum(wps2/(w*(w-wcs)))
    # P = 1-sum(wps2/w**2)
    S = 1/2*(R+L)
    D = 1/2*(R-L)

    return S,D,P,R,L

def suprathermal_distribution_bell(vperp,vpar):
    a=4.9e5
    b=8.3e14
    c=5.4e23

    v0=1

    v=100*np.sqrt(vperp*vperp+vpar*vpar+v0*v0)  #cm/s
    # vm=np.sqrt(vperp*vperp+vpar*vpar+v0*v0)  #m/s

    # v=100*(vm/(np.sqrt(1-(vm**2/const.c_light**2))))

    fbell=(a/(v**4))-(b/(v**5))+(c/(v**6))


    ####
    #Bortnik, J., Inan, U. S., and Bell, T. F. (2006), Landau damping and resultant unidirectional 
    # propagation of chorus waves, Geophys. Res. Lett., 33, L03102, doi:10.1029/2005GL024553.  --> 
    # 
    # f=10*fbell
    #####
    f=fbell*(10**12)  #s^3/m^6 from s^3/cm^6

    return f

def spatialdamping(f, kperp, kpar, w, m, wch, qh, mh,qs, Ns, ms, nus, B0):

# Compute the spatial damping rate ki in a hot plasma

    theta=np.arctan2(kperp,kpar)
    n=np.sqrt((const.c_light**2/w**2)*(kperp**2+kpar**2))

    [S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0)

    

    A = S*np.sin(theta)**2 + P*np.cos(theta)**2
    B = R*L*np.sin(theta)**2 + P*S*(1+np.cos(theta)**2)

    # print('A=',A)
    # print('B=',B)
    ki = 0
 
    for ii in range(0,len(f)):
        Di=hot_dispersion_imag(f[ii], kperp, kpar, w, m , wch[ii],
                                    qh[ii], mh[ii], qs , Ns, ms, nus, B0)
        ki=ki + -(w/const.c_light)*(1/2)*(1/(4*n*(2*A*(n**2)-B)))*Di
        # ki=ki + -(w/const.c_light)*(1/(4*n*(2*A*(n**2)-B)))*Di #ST removed (1/2) from ki

    return ki

def integrand(f,vperp,kperp,kpar,w,m,wch,qh,mh,R,L,P,S):

    theta = np.arctan2(kperp,kpar)
    #refractive index
    n=np.sqrt((const.c_light**2/w**2)*(kperp**2+kpar**2))

    ct=np.cos(theta)
    st=np.sin(theta)

    vperp_vec=vperp
    ret=np.zeros(np.size(vperp_vec))

    


    # for vperp_i in range (0,len(vperp_vec)):

    for vperp_i in range (0,1):
        # vperp=vperp_vec[vperp_i]
        vperp=vperp_vec

        suminteg=0

        for mi in range(0,len(m)):
            Jm=scp.jv(m[mi], kperp*vperp/wch)
            Jm_m1=scp.jv(m[mi]-1, kperp*vperp/wch)
            Jm_p1=scp.jv(m[mi]+1, kperp*vperp/wch)
            

            Jm=np.real(Jm)
            Jm_m1=np.real(Jm_m1)
            Jm_p1=np.real(Jm_p1)

            #evaluate vpar only at resonance
            vpar = (w-m[mi]*wch)/kpar
            # evaluate G1
            G1 = fG1(f, vperp, vpar, kperp, kpar, w)
            # evaluate G2
            G2 = fG2(f, vperp, vpar, kperp, kpar, w, m[mi], wch)      

            suminteg=suminteg+  \
                (G1*((P-n**2*st**2)*(2*(L-n**2)*vperp*Jm_p1*Jm_p1+2*vperp*(R-n**2)*Jm_m1**2+    \
                    n**2*st*st*vperp*(Jm_p1-Jm_m1)*(Jm_p1-Jm_m1))-    \
                    n**2*ct*st*(2*vpar*Jm*(Jm_p1*(R-n**2)+Jm_m1*(L-n**2))+  \
                    n**2*ct*st*vperp*(Jm_p1-Jm_m1)*(Jm_p1-Jm_m1)))+  \
                    G2*(4*vpar*Jm*((L-n**2)*(R-n**2)+n*n*st*st*(S-n**2))-
                    2*n*n*ct*st*((R-n**2)*vperp*Jm_m1+(L-n**2)*vperp*Jm_p1)))
            
        ret[vperp_i]=-2*np.pi*np.pi*((qh*qh/mh/const.eps0)/(w*np.absolute(kpar)))*suminteg*vperp

    return ret
def fG1(f, vperp, vpar, kperp, kpar, w):

    # f = distribution function, a function handle that takes as input 
    #     (vperp,vpar) and returns the value of the distribution at that
    #     point.
    # vperp:  perpendicular velocity
    # vpar:   parallel velocity
    # kperp:  perpendicular component of k vector
    # kpar:   parallel component of k vector
    # w:      frequency, rad/s 

    eps=2.220446049250313e-16
    DEL=1e-8
    d=DEL*np.absolute(vperp)
    if d<10*eps:
        d=10*eps 

    df_dvperp = (f(vperp+d, vpar)-f(vperp-d,vpar))/(2*d)

    d=DEL*np.absolute(vpar)
    if d<10*eps:
        d=10*eps 

    df_dvpar = (f(vperp, vpar+d)-f(vperp,vpar-d))/(2*d)

    G1 = df_dvperp-(kpar/w)*(vpar*df_dvperp - vperp*df_dvpar)

    return G1

def fG2(f, vperp, vpar, kperp, kpar, w, m, wch):

    # f = distribution function, a function handle that takes as input 
    #     (vperp,vpar) and returns the value of the distribution at that
    #     point.
    # vperp:  perpendicular velocity
    # vpar:   parallel velocity
    # kperp:  perpendicular component of k vector
    # kpar:   parallel component of k vector
    # w:      frequency, rad/s 
    # m:      resonance (0=landau, -1,+1=cyclotron, etc.)
    # wch:    hot (signed) gyrofrequency, negative for electrons
    DEL=1e-8

    eps=2.220446049250313e-16
    d=DEL*np.absolute(vperp)
    if d<10*eps:
        d=10*eps 

    df_dvperp = (f(vperp+d, vpar)-f(vperp-d,vpar))/(2*d)

    d=DEL*np.absolute(vpar)
    if d<10*eps:
        d=10*eps 

    df_dvpar = (f(vperp, vpar+d)-f(vperp,vpar-d))/(2*d)

    Jm=np.real(scp.jv(m,kperp*vperp/wch))

    # Should this be -(m*wch) or +(m*wch)???
    # G2 = Jm*(df_dvpar-(m*wch+eps)/(w*vperp+eps)*(vpar*df_dvperp - vperp*df_dvpar))
    G2 = Jm*(df_dvpar-(m*wch+eps)/(w*vperp+eps)*(vpar*df_dvperp - vperp*df_dvpar))
    return G2

def hot_dispersion_imag(f, kperp, kpar, w, m, wch, qh, mh, qs, Ns, ms, nus, B0 ):

    eps=2.220446049250313e-16

    [S,D,P,R,L] = stix_parameters(w, qs, Ns, ms, nus, B0)
    

    integrand_vperp=lambda vperp: integrand(f, vperp, kperp, kpar, w, m,
                                      wch, qh, mh,R, L, P, S)
    

    k=const.c_light
    integrand_vperpnorm= lambda vperp: k*integrand_vperp(k*vperp)

    integrand_t=lambda t: ((1+eps)/(t**2+eps))*integrand_vperpnorm((1-t+eps)/(t+eps))

    # integrated_integrand=scint.quad(integrand_t,0,1,epsabs=1e-3,maxp1=100,limit=100)
    integrated_integrand=scint.quad(integrand_vperpnorm,0,1,epsabs=1e-3,maxp1=100,limit=100)
    # integrated_integrand=scint.quad(integrand_vperp,0,k,maxp1=100,limit=100)
    ret=integrated_integrand[0]
  
    return ret




def read_input_ray(ray_file_name):


    ray_in=pd.read_csv(ray_file_name,delim_whitespace=True, header=None)
    ray_in.columns=['ray','ray_stop','time','posx','posy','posz','vprelx','vprely','vprelz',
                    'vgrelx','vgrely','vgrelz','nx','ny','nz','Bx','By','Bz','w','Nspec','qs1','qs2', 'qs3',
                    'qs4', 'ms1','ms2','ms3','ms4','Ns1','Ns2','Ns3','Ns4','nus1','nus2','nus3','nus4' ]

    timef=ray_in.time
    posxf=ray_in.posx
    posyf=ray_in.posy
    poszf=ray_in.posz
    vprelxf=ray_in.vprelx
    vprelyf=ray_in.vprely
    vprelzf=ray_in.vprelz
    vgrelxf=ray_in.vgrelx
    vgrelyf=ray_in.vgrely
    vgrelzf=ray_in.vgrelz
    nxf=ray_in.nx
    nyf=ray_in.ny
    nzf=ray_in.nz
    Bxf=ray_in.Bx
    Byf=ray_in.By
    Bzf=ray_in.Bz
    wf=ray_in.w
    Nspecf=ray_in.Nspec
    qs1f=ray_in.qs1
    qs2f=ray_in.qs2
    qs3f=ray_in.qs3
    qs4f=ray_in.qs4
    ms1f=ray_in.ms1
    ms2f=ray_in.ms2
    ms3f=ray_in.ms3
    ms4f=ray_in.ms4
    Ns1f=ray_in.Ns1
    Ns2f=ray_in.Ns2
    Ns3f=ray_in.Ns3
    Ns4f=ray_in.Ns4
    nus1f=ray_in.nus1
    nus2f=ray_in.nus2
    nus3f=ray_in.nus3
    nus4f=ray_in.nus4

    freq=wf/(2*np.pi)

    data={'time':timef,'posx':posxf,'posy':posyf,'posz':poszf,'vprelx':vprelxf,'vprely':vprelyf,'vprelz':vprelzf,
                   'vgrelx':vgrelxf,'vgrely':vgrelyf,'vgrelz':vgrelzf,'nx':nxf,'ny':nyf,'nz':nzf,'Bx':Bxf,'By':Byf,'Bz':Bzf,'w':wf,
                   'Nspec':Nspecf,'qs1':qs1f,'qs2':qs2f, 'qs3':qs3f,
                   'qs4':qs4f, 'ms1':ms1f,'ms2':ms2f,'ms3':ms3f,'ms4':ms4f,'Ns1':Ns1f,'Ns2':Ns2f,'Ns3':Ns3f,'Ns4':Ns4f,'nus1':nus1f,
                   'nus2':nus2f,'nus3':nus3f,'nus4':nus4f,'freq':freq}    
    df = pd.DataFrame(data)

    return df

def landau_damping(ray_file):

    ray_in=read_input_ray(ray_file)

    # print(df.posx)
    psi=np.zeros(len(ray_in.time))
    freq=ray_in.freq

    for kk in range (0,len(ray_in.time)):

        nmag0=np.sqrt(ray_in.nx[kk]*ray_in.nx[kk]+ray_in.ny[kk]*ray_in.ny[kk]+ray_in.nz[kk]*ray_in.nz[kk])
        Bmag0=np.sqrt(ray_in.Bx[kk]*ray_in.Bx[kk]+ray_in.By[kk]*ray_in.By[kk]+ray_in.Bz[kk]*ray_in.Bz[kk])
        pmag0=np.sqrt(ray_in.posx[kk]*ray_in.posx[kk]+ray_in.posy[kk]*ray_in.posy[kk]+ray_in.posz[kk]*ray_in.posz[kk])
        n=[ray_in.nx[kk]/nmag0,ray_in.ny[kk]/nmag0,ray_in.nz[kk]/nmag0]
        B=[ray_in.Bx[kk]/Bmag0,ray_in.By[kk]/Bmag0,ray_in.Bz[kk]/Bmag0]
        p=[ray_in.posx[kk]/pmag0,ray_in.posy[kk]/pmag0,ray_in.posz[kk]/pmag0]

        nBdot=np.dot(n,B)
        ndot=np.dot(n,n)
        Bdot=np.dot(B,B)

        cos2psi=(nBdot*nBdot)/(ndot*Bdot)
        cospsi=np.sqrt(cos2psi)
        psi2=np.arccos(cospsi)
        psi[kk]=psi2
        # print(np.rad2deg(psi[kk]))
        
    theta_res=np.zeros(len(ray_in.time))
    Y=np.zeros(len(ray_in.time))
    S_stix=np.zeros(len(ray_in.time))
    D_stix=np.zeros(len(ray_in.time))
    P_stix=np.zeros(len(ray_in.time))
    R_stix=np.zeros(len(ray_in.time))
    L_stix=np.zeros(len(ray_in.time))
    wce=np.zeros(len(ray_in.time))
    for ll in range(0,len(ray_in.time)):
        Ns=[ray_in.Ns1[ll],ray_in.Ns2[ll],ray_in.Ns3[ll],ray_in.Ns4[ll]]
        qs=[ray_in.qs1[ll],ray_in.qs2[ll],ray_in.qs3[ll],ray_in.qs4[ll]]
        ms=[ray_in.ms1[ll],ray_in.ms2[ll],ray_in.ms3[ll],ray_in.ms4[ll]]
        B0mag=np.sqrt(ray_in.Bx[ll]*ray_in.Bx[ll]+ray_in.By[ll]*ray_in.By[ll]+ray_in.Bz[ll]*ray_in.Bz[ll])
        w=ray_in.w[ll]



        # print(Ns[0],qs[0],ms[0],w)
        wps21 = (Ns[0]*qs[0]**2/ms[0]/const.eps0)
        wcs1 = ((qs[0]*B0mag)/ms[0])

        wps22 = (Ns[1]*qs[1]**2/ms[1]/const.eps0)
        wcs2 = ((qs[1]*B0mag)/ms[1])

        wps23 = (Ns[2]*qs[2]**2/ms[2]/const.eps0)
        wcs3 = ((qs[2]*B0mag)/ms[2])

        wps24 = (Ns[3]*qs[3]**2/ms[3]/const.eps0)
        wcs4 = ((qs[3]*B0mag)/ms[3])

        # S_,D_,P_,R_,L_=stix_parameters(w, ray_in.Ns1[ll], ray_in.Ns2[ll], ray_in.Ns3[ll], ray_in.Ns4[ll], B0mag)

        # print(ms[0],const.me)

        R_=1-(wps21/(w*(w+wcs1)))-(wps22/(w*(w+wcs2)))-(wps23/(w*(w+wcs3)))-(wps24/(w*(w+wcs4)))
        L_=1-(wps21/(w*(w-wcs1)))-(wps22/(w*(w-wcs2)))-(wps23/(w*(w-wcs3)))-(wps24/(w*(w-wcs4)))
        P_=1-(wps21/w**2)-(wps22/w**2)-(wps23/w**2)-(wps24/w**2)
        S_ = 1/2*(R_+L_)
        D_ = 1/2*(R_-L_)


        # thetares=res_angle(P_,S_)
        # thetares3=np.pi-np.arccos(w/wcs1)
        tantheta_res2=-P_/S_
        tan_theta=np.sqrt(np.absolute(tantheta_res2))
        theta_res_tmp=np.arctan(tan_theta)
        # print(P_/S_,np.absolute(P_)/np.absolute(S_))

        theta_res_deg_tmp=np.rad2deg(theta_res_tmp)
        theta_res[ll]=theta_res_deg_tmp

        Ytmp=np.absolute((np.rad2deg(psi[ll])-theta_res_deg_tmp))
        # print(Ytmp)
        Y[ll]=Ytmp
        S_stix[ll]=S_
        D_stix[ll]=D_
        P_stix[ll]=P_
        R_stix[ll]=R_
        L_stix[ll]=L_
        wce[ll]=wcs1

    raytime = dt.datetime(2010,1,1,0,0,0)
    tmp_coords = coord.Coords(list(zip(ray_in.posx, ray_in.posy, ray_in.posz)),'SM','car',units=['m','m','m'])
    tvec_datetime = [raytime + dt.timedelta(seconds=s) for s in ray_in.time]
    tmp_coords.ticks = Ticktock(tvec_datetime)
    tmp_coords = tmp_coords.convert('MAG','sph')

    alt_ray=tmp_coords.radi
    lat_ray=tmp_coords.lati
    lon_ray=tmp_coords.long

    L_r = (alt_ray/const.Re)/(np.cos(np.deg2rad(lat_ray))*np.cos(np.deg2rad(lat_ray)))

    ref_ind=np.sqrt(ray_in.nx**2+ray_in.ny**2+ray_in.nz**2)
    eta_sq_plus,eta_sq_minus=ref_index(psi,S_stix,P_stix,R_stix,L_stix)



    t=ray_in.time
    w=ray_in.w[kk]
    # nus=[ray_in.nus1,ray_in.nus2,ray_in.nus3]
    nus=[0,0,0]



    fe= lambda vperp,vpar: suprathermal_distribution_bell(vperp,vpar)


    #resonance orders to include in damping (m=0 Landau is the most signifficant by far)
    mres=[-2,-1,0,1,2]

    fs = [fe]    

    magnitude=np.zeros(np.size(t))
    kis=np.zeros(np.size(t))
    magnitude[0]=1  #set initial magnitude of the power, normalised thus magnitude[0]=1
            
    print("step of total steps:")
    for ii in range(1,len(t)):
        print('%d of %d' %(ii-1, len(t)-1), end='\r')
        # qe_h = const.qe
        # me_h=const.me
        normB=np.sqrt(ray_in.Bx[ii]*ray_in.Bx[ii]+ray_in.By[ii]*ray_in.By[ii]+ray_in.Bz[ii]*ray_in.Bz[ii])
        wce_h=((const.qe*normB)/const.me)
        # vector of hot plasma properties (one per hot species)
        wchs=[wce_h]
        qhs=[const.qe]
        mhs=[const.me]

        kx=ray_in.nx[ii]*w/const.c_light
        ky=ray_in.ny[ii]*w/const.c_light
        kz=ray_in.nz[ii]*w/const.c_light
        k=np.array([kx,ky,kz])
        kmag=np.sqrt(kx*kx+ky*ky+kz*kz)

        Bhat=np.array([ray_in.Bx[ii]/normB,ray_in.By[ii]/normB,ray_in.Bz[ii]/normB])


        vgrel=[ray_in.vgrelx[ii],ray_in.vgrely[ii],ray_in.vgrelz[ii]]
        kpar=np.dot(k,Bhat)

        kperp=k-kpar*Bhat

        kperp=np.sqrt(kperp[0]*kperp[0]+kperp[1]*kperp[1]+kperp[2]*kperp[2])
        kmag=np.sqrt(kx*kx+ky*ky+kz*kz)
        vgrelmag=np.sqrt(ray_in.vgrelx[ii]*ray_in.vgrelx[ii]+ray_in.vgrely[ii]*ray_in.vgrely[ii]+ray_in.vgrelz[ii]*ray_in.vgrelz[ii])

        qs=np.array([ray_in.qs1[ii],ray_in.qs2[ii],ray_in.qs3[ii],ray_in.qs4[ii]])
        Ns=np.array([ray_in.Ns1[ii],ray_in.Ns2[ii],ray_in.Ns3[ii],ray_in.Ns4[ii]])
        ms=np.array([ray_in.ms1[ii],ray_in.ms2[ii],ray_in.ms3[ii],ray_in.ms4[ii]])
        nus=np.array([ray_in.nus1[ii],ray_in.nus2[ii],ray_in.nus3[ii],ray_in.nus4[ii]])


        if kmag!=0:
            ki=spatialdamping(fs,kperp,kpar,w,mres,wchs,qhs,mhs,
                            qs,Ns,ms,nus, normB)
            ki_along_vg=ki*((np.dot(k,vgrel))/(kmag*vgrelmag))

            posx=ray_in.posx[ii]-ray_in.posx[ii-1]
            posy=ray_in.posy[ii]-ray_in.posy[ii-1]
            posz=ray_in.posz[ii]-ray_in.posz[ii-1]

            dist=np.sqrt(posx**2+posy**2+posz**2)
            # print('dist',dist)
            kis[ii]=ki_along_vg
            # magnitude[ii]=magnitude[ii-1]*np.exp(-2*dist*ki_along_vg)
            magnitude[ii]=magnitude[ii-1]*np.exp(-dist*ki) #ST changed ki_along_vg to ki
            if magnitude[ii]<0.01:
                break
            # print('ki=',ki)
            # print('kvg=',ki_along_vg)
            # print('kiss=',k)
    #             print('mag=',magnitude[ii])
        else:
            print('Re[n]=0, not solving evanescent mode')
    print(kis[1],magnitude[0],magnitude[2])


    print('done')
    print ("Finished with Landau damping")


    mag_zero=np.where(magnitude == 0)

    print(mag_zero[0][0])

    # --------------- Latex Plot Beautification --------------------------
    print('Start plotting...')
    fig_width_pt = 650.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = 6  # width in inches
    fig_height = 6      # height in inches
    fig_size =  [fig_width+1,fig_height+1]
    params = {'backend': 'ps',
                'axes.labelsize': 14,
                'font.size': 14,
                'legend.fontsize': 14,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'text.usetex': False,
                'figure.figsize': fig_size}
    plt.rcParams.update(params)

    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)
    #--------------- Latex Plot Beautification --------------------------
    # plt.style.use('ggplot')

    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    ax2.set_ylim(0,1)
    ax2.set_xlim(0,t[mag_zero[0][0]])

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Normalised Wave Power')

    ax2.plot(t[0:mag_zero[0][0]],magnitude[0:mag_zero[0][0]],c=tableau20[4], alpha=0.75)
    plt.grid()
    plt.show()

    print('Saving file...')
    psi_file=np.round(np.rad2deg(psi[0]),2) 
    freq_file=np.round(freq[0],0)     
    L_file=np.round(L_r[0],2) 
    lat_file=np.round(lat_ray[0],2)


    data={'time':ray_in.time[0:mag_zero[0][0]],
            'posx':ray_in.posx[0:mag_zero[0][0]],'posy':ray_in.posy[0:mag_zero[0][0]],'posz':ray_in.posz[0:mag_zero[0][0]],
            'vprelx':ray_in.vprelx[0:mag_zero[0][0]],'vprely':ray_in.vprely[0:mag_zero[0][0]],'vprelz':ray_in.vprelz[0:mag_zero[0][0]],
            'vgrelx':ray_in.vgrelx[0:mag_zero[0][0]],'vgrely':ray_in.vgrely[0:mag_zero[0][0]],'vgrelz':ray_in.vgrelz[0:mag_zero[0][0]],
            'nx':ray_in.nx[0:mag_zero[0][0]],'ny':ray_in.ny[0:mag_zero[0][0]],'nz':ray_in.nz[0:mag_zero[0][0]],
            'Bx':ray_in.Bx[0:mag_zero[0][0]],'By':ray_in.By[0:mag_zero[0][0]],'Bz':ray_in.Bz[0:mag_zero[0][0]],
            'w':ray_in.w[0:mag_zero[0][0]],
            'qs1':ray_in.qs1[0:mag_zero[0][0]],'qs2':ray_in.qs2[0:mag_zero[0][0]], 'qs3':ray_in.qs3[0:mag_zero[0][0]],'qs4':ray_in.qs4[0:mag_zero[0][0]],
            'ms1':ray_in.ms1[0:mag_zero[0][0]],'ms2':ray_in.ms2[0:mag_zero[0][0]],'ms3':ray_in.ms3[0:mag_zero[0][0]],'ms4':ray_in.ms4[0:mag_zero[0][0]],
            'Ns1':ray_in.Ns1[0:mag_zero[0][0]],'Ns2':ray_in.Ns2[0:mag_zero[0][0]],'Ns3':ray_in.Ns3[0:mag_zero[0][0]],'Ns4':ray_in.Ns4[0:mag_zero[0][0]],
            'psi':np.rad2deg(psi[0:mag_zero[0][0]]),'theta_res':theta_res[0:mag_zero[0][0]],'Y':Y[0:mag_zero[0][0]],
            'L':L_r[0:mag_zero[0][0]],'alt':alt_ray[0:mag_zero[0][0]],'lat':lat_ray[0:mag_zero[0][0]],
            'lon':lon_ray[0:mag_zero[0][0]],'damp':magnitude[0:mag_zero[0][0]],
            'S_stix':S_stix[0:mag_zero[0][0]],'D_stix':D_stix[0:mag_zero[0][0]],'P_stix':P_stix[0:mag_zero[0][0]],
            'R_stix':R_stix[0:mag_zero[0][0]],'L_stix':L_stix[0:mag_zero[0][0]],'wce':wce[0:mag_zero[0][0]]}    
    df = pd.DataFrame(data)


    path_csv='example_rays/Landau_output/'
    # Create directory for inputs/outputs if doesn't already exist
    if not os.path.exists(path_csv):
        os.makedirs(path_csv)

    df.to_csv(path_csv+'L%d_freq%d_psi%d_lat_%d_damping.csv' %(L_file,freq_file,psi_file,lat_file),sep=',',index=False)
    
    print('L%d_freq%d_psi%d_lat_%d_damping.csv' %(L_file,freq_file,psi_file,lat_file))