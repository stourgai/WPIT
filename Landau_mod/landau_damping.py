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
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import environment_mod as env
import waveproperties_mod as wave
import ray_mod as ray
import Landau_mod as landau

def landau_damping(ray_file,distr):

    ray_in=ray.read_input_ray(ray_file)

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
        fe= lambda vperp,vpar: landau.bell_distribution(vperp,vpar)

    elif distr=='Bortnik':
        fe= lambda vperp,vpar: landau.bortnik_distribution(vperp,vpar)
    
    elif distr=='Bimaxw':
        fe= lambda vperp,vpar: landau.bi_maxwellian_distribution(vperp,vpar)

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
            ki=landau.spatialdamping(fs,kperp,kpar,w,mres,wchs,qhs,mhs,
                            qs,Ns,ms,nus, normB)
            ki_along_vg=ki*((np.dot(k,vgrel))/(kmag*vgrelmag))

            posx=ray_in.posx[ii]-ray_in.posx[ii-1]
            posy=ray_in.posy[ii]-ray_in.posy[ii-1]
            posz=ray_in.posz[ii]-ray_in.posz[ii-1]

            dist=np.sqrt(posx**2+posy**2+posz**2)
            kis[ii]=ki_along_vg
            magnitude[ii]=magnitude[ii-1]*np.exp(-dist*ki) #ST changed ki_along_vg to ki
            if magnitude[ii]<0.01:
                break
        else:
            print('Re[n]=0, not solving evanescent mode')
    print(kis[1],magnitude[0],magnitude[2])


    print('done')
    print ("Finished with Landau damping")


    mag_zero=np.where(magnitude == 0)

    print(mag_zero[0][0])

    
    fig2, ax2 = plt.subplots(figsize=(9,10))
    # ax2=fig2.add_subplot(111)

    ax2.set_ylim(0,1)
    ax2.set_xlim(0,t[mag_zero[0][0]])

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Normalised Wave Power')

    ax2.plot(t[0:mag_zero[0][0]],magnitude[0:mag_zero[0][0]],c='tab:green', alpha=0.75)
    plt.grid()
    plt.show()

    print('Saving file...')



    data={'time':ray_in.time[0:mag_zero[0][0]],'damp':magnitude[0:mag_zero[0][0]]}    
    df = pd.DataFrame(data)


    path_csv='example_rays/Landau_output/'
    # Create directory for inputs/outputs if doesn't already exist
    if not os.path.exists(path_csv):
        os.makedirs(path_csv)

    df.to_csv(path_csv+'damping.csv',sep=',',index=False)
    
 