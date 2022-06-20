import pandas as pd    
import numpy as np
import matplotlib.pyplot as plt
import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
from WPIT.LandauDamp_mod.RayUtils_mod import read_input_ray
from WPIT.LandauDamp_mod import distribution_bell,distribution_bortnik,distribution_bimaxwellian,distribution_golden2
from scipy import integrate as scint
from scipy import special as scp


def spatialdamping(f, kperp, kpar, w, m, wch, qh, mh,qs, Ns, ms, nus, B0):

    # Compute the spatial damping rate ki in a hot plasma

    theta=np.arctan2(kperp,kpar)
    n=np.sqrt((env.const.c_light**2/w**2)*(kperp**2+kpar**2))

    [S,D,P,R,L] = wave.stix_parameters(w, Ns[0],Ns[1],Ns[2],Ns[3], B0)

    

    A = S*np.sin(theta)**2 + P*np.cos(theta)**2
    B = R*L*np.sin(theta)**2 + P*S*(1+np.cos(theta)**2)

    # print('A=',A)
    # print('B=',B)
    ki = 0
 
    for ii in range(0,len(f)):
        Di=hot_dispersion_imag(f[ii], kperp, kpar, w, m , wch[ii],
                                    qh[ii], mh[ii], qs , Ns, ms, nus, B0)
        ki=ki + -(w/env.const.c_light)*(1/2)*(1/(4*n*(2*A*(n**2)-B)))*Di
        # ki=ki + -(w/const.c_light)*(1/(4*n*(2*A*(n**2)-B)))*Di #ST removed (1/2) from ki

    return ki

def integrand(f,vperp,kperp,kpar,w,m,wch,qh,mh,R,L,P,S):

    theta = np.arctan2(kperp,kpar)
    #refractive index
    n=np.sqrt((env.const.c_light**2/w**2)*(kperp**2+kpar**2))

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
            
        ret[vperp_i]=-2*np.pi*np.pi*((qh*qh/mh/env.const.eps0)/(w*np.absolute(kpar)))*suminteg*vperp

    return ret

def hot_dispersion_imag(f, kperp, kpar, w, m, wch, qh, mh, qs, Ns, ms, nus, B0 ):

    eps=2.220446049250313e-16

    [S,D,P,R,L] = wave.stix_parameters(w, Ns[0],Ns[1],Ns[2],Ns[3], B0)
    

    integrand_vperp=lambda vperp: integrand(f, vperp, kperp, kpar, w, m,
                                      wch, qh, mh,R, L, P, S)
    

    k=env.const.c_light
    integrand_vperpnorm= lambda vperp: k*integrand_vperp(k*vperp)

    integrand_t=lambda t: ((1+eps)/(t**2+eps))*integrand_vperpnorm((1-t+eps)/(t+eps))

    # integrated_integrand=scint.quad(integrand_t,0,1,epsabs=1e-3,maxp1=100,limit=100)
    integrated_integrand=scint.quad(integrand_vperpnorm,0,1,epsabs=1e-3,maxp1=100,limit=100)
    # integrated_integrand=scint.quad(integrand_vperp,0,k,maxp1=100,limit=100)
    ret=integrated_integrand[0]
  
    return ret



def fG1(f, vperp, vpar, kperp, kpar, w):

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

    G2 = Jm*(df_dvpar-(m*wch+eps)/(w*vperp+eps)*(vpar*df_dvperp - vperp*df_dvpar))
    return G2

def landau_damping(ray_file,distr):

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
        


    t=ray_in.time
    w=ray_in.w[kk]
    # nus=[ray_in.nus1,ray_in.nus2,ray_in.nus3]
    nus=[0,0,0]


    if distr=='Bell':
        fe= lambda vperp,vpar: distribution_bell(vperp,vpar)

    elif distr=='Bortnik':
        fe= lambda vperp,vpar: distribution_bortnik(vperp,vpar)
    
    elif distr=='Bimaxw':
        fe= lambda vperp,vpar: distribution_bimaxwellian(vperp,vpar)

    elif distr=='Golden':
        fe= lambda vperp,vpar: distribution_golden2(vperp,vpar)

    else:
        print('No ',distr,' distribution. Avalaible options: Bell, Bortnik or Bimaxw')


    #resonance orders to include in damping (m=0 Landau is the most signifficant by far)
    mres=[-2,-1,0,1,2]

    fs = [fe]    

    magnitude=np.zeros(np.size(t))
    kis=np.zeros(np.size(t))
    magnitude[0]=1  #set initial magnitude of the power, normalised thus magnitude[0]=1
            
    print("step of total steps:")
    for ii in range(1,len(t)):
        print('%d of %d' %(ii-1, len(t)-1), end='\r')
        # qe_h = env.const.qe
        # me_h=env.const.me
        normB=np.sqrt(ray_in.Bx[ii]*ray_in.Bx[ii]+ray_in.By[ii]*ray_in.By[ii]+ray_in.Bz[ii]*ray_in.Bz[ii])
        wce_h=((env.const.qe*normB)/env.const.me)
        # vector of hot plasma properties (one per hot species)
        wchs=[wce_h]
        qhs=[env.const.qe]
        mhs=[env.const.me]

        kx=ray_in.nx[ii]*w/env.const.c_light
        ky=ray_in.ny[ii]*w/env.const.c_light
        kz=ray_in.nz[ii]*w/env.const.c_light
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
            ki=spatialdamping(fs,kperp,kpar,w,mres,wchs,qhs,mhs,qs,Ns,ms,nus, normB)
            ki_along_vg=ki*((np.dot(k,vgrel))/(kmag*vgrelmag))

            posx=ray_in.posx[ii]-ray_in.posx[ii-1]
            posy=ray_in.posy[ii]-ray_in.posy[ii-1]
            posz=ray_in.posz[ii]-ray_in.posz[ii-1]

            dist=np.sqrt(posx**2+posy**2+posz**2)
            kis[ii]=ki_along_vg
            magnitude[ii]=magnitude[ii-1]*np.exp(-dist*ki) 
            # if magnitude[ii]<0.01:
            #     break
        else:
            print('Re[n]=0, not solving evanescent mode')
    print(kis[1],magnitude[0],magnitude[2])


    print('done')
    print ("Finished with Landau damping")


    # mag_zero=np.where(magnitude == 0)

    # print(mag_zero[0][0])

    
    fig2, ax2 = plt.subplots(figsize=(9,10))
    # ax2=fig2.add_subplot(111)

    ax2.set_ylim(0,1)
    # ax2.set_xlim(0,t[mag_zero[0][0]])

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Normalised Wave Power')

    ax2.plot(t,magnitude,c='tab:green', alpha=0.75)
    plt.grid()
    plt.show()

    print('Saving file...')



    data={'time':ray_in.time,'damp':magnitude}    
    df = pd.DataFrame(data)



    
    df.to_csv(ray_file+'_damping.csv',sep=',',index=False)

    print('Saved file as:', ray_file+'_damping.csv')