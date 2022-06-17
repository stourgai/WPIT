import pandas as pd    
import numpy as np
import datetime as dt
from spacepy import coordinates as coord
from spacepy.time import Ticktock

import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
from .read_input_ray import read_input_ray

def append_ray(ray_file_name):

    ray_=read_input_ray(ray_file_name)
    damp_file=pd.read_csv(ray_file_name+'_damping.csv')
    damping=damp_file['damp']

    damp_cut=np.zeros(len(damping))

    for i in range(0,len(damping)):
        if damping[i]>0.01:
            damp_cut[i]=damping[i]

        else:
            break


    time=ray_.time[0:len(damp_cut)]

    posx=ray_.posx[0:len(damp_cut)]
    posy=ray_.posy[0:len(damp_cut)]
    posz=ray_.posz[0:len(damp_cut)]
    vprelx=ray_.vprelx[0:len(damp_cut)]
    vprely=ray_.vprely[0:len(damp_cut)]
    vprelz=ray_.vprelz[0:len(damp_cut)]
    vgrelx=ray_.vgrelx[0:len(damp_cut)]
    vgrely=ray_.vgrely[0:len(damp_cut)]
    vgrelz=ray_.vgrelz[0:len(damp_cut)]
    nx=ray_.nx[0:len(damp_cut)]
    ny=ray_.ny[0:len(damp_cut)]
    nz=ray_.nz[0:len(damp_cut)]
    Bx=ray_.Bx[0:len(damp_cut)]
    By=ray_.By[0:len(damp_cut)]
    Bz=ray_.Bz[0:len(damp_cut)]
    qs1=ray_.qs1[0:len(damp_cut)]
    qs2=ray_.qs2[0:len(damp_cut)]
    qs3=ray_.qs3[0:len(damp_cut)]    
    qs4=ray_.qs4[0:len(damp_cut)]
    ms1=ray_.ms1[0:len(damp_cut)]
    ms2=ray_.ms2[0:len(damp_cut)]
    ms3=ray_.ms3[0:len(damp_cut)]    
    ms4=ray_.ms4[0:len(damp_cut)]
    Ns1=ray_.Ns1[0:len(damp_cut)]
    Ns2=ray_.Ns2[0:len(damp_cut)]
    Ns3=ray_.Ns3[0:len(damp_cut)]    
    Ns4=ray_.Ns4[0:len(damp_cut)]
    w=ray_.w[0:len(damp_cut)]
    freq=w/(2*np.pi)

    psi=np.zeros(len(time))
    theta_res=np.zeros(len(time))
    Y=np.zeros(len(time))
    S_stix=np.zeros(len(time))
    D_stix=np.zeros(len(time))
    P_stix=np.zeros(len(time))
    R_stix=np.zeros(len(time))
    L_stix=np.zeros(len(time))

    wce_out=np.zeros(len(time))
    wcH_out=np.zeros(len(time))
    wcHe_out=np.zeros(len(time))
    wcO_out=np.zeros(len(time))
    wpe_out=np.zeros(len(time))
    wpH_out=np.zeros(len(time))
    wpHe_out=np.zeros(len(time))
    wpO_out=np.zeros(len(time))    
    wlhr_out=np.zeros(len(time))
    wuhr_out=np.zeros(len(time))
    ref_out=np.zeros(len(time))
    gendrin_out=np.zeros(len(time))

    for kk in range (0,len(time)):

        nmag0=np.sqrt(nx[kk]*nx[kk]+ny[kk]*ny[kk]+nz[kk]*nz[kk])
        Bmag0=np.sqrt(Bx[kk]*Bx[kk]+By[kk]*By[kk]+Bz[kk]*Bz[kk])
        pmag0=np.sqrt(posx[kk]*posx[kk]+posy[kk]*posy[kk]+posz[kk]*posz[kk])
        nhat=[nx[kk]/nmag0,ny[kk]/nmag0,nz[kk]/nmag0]
        Bhat=[Bx[kk]/Bmag0,By[kk]/Bmag0,Bz[kk]/Bmag0]
        phat=[posx[kk]/pmag0,posy[kk]/pmag0,posz[kk]/pmag0]

        nBdot=np.dot(nhat,Bhat)
        ndot=np.dot(nhat,nhat)
        Bdot=np.dot(Bhat,Bhat)

        cos2psi=(nBdot*nBdot)/(ndot*Bdot)
        cospsi=np.sqrt(cos2psi)
        psitmp=np.arccos(cospsi)

        wce=env.omega_cyclotron(Bmag0,env.const.qe,env.const.me)
        wcH=env.omega_cyclotron(Bmag0,env.const.qi,env.const.mH)
        wcHe=env.omega_cyclotron(Bmag0,env.const.qi,env.const.mHe)
        wcO=env.omega_cyclotron(Bmag0,env.const.qi,env.const.mO)

        wpe=env.omega_plasma(Ns1[kk],env.const.qe,env.const.me)
        wpH=env.omega_plasma(Ns2[kk],env.const.qi,env.const.mH)
        wpHe=env.omega_plasma(Ns3[kk],env.const.qe,env.const.mHe)
        wpO=env.omega_plasma(Ns4[kk],env.const.qi,env.const.mO)

        wlhr=env.omega_lhr(wce,wpe,wcH,wpH)

        wuhr=env.omega_uhr(wce,wpe)

        S_stixtmp,D_stixtmp,P_stixtmp,R_stixtmp,L_stixtmp=wave.stix_parameters(w[kk], Ns1[kk], Ns2[kk], Ns3[kk], Ns4[kk], Bmag0)
        thetarestmp=wave.res_angle(P_stixtmp,S_stixtmp)
        theta_res_deg=np.rad2deg(thetarestmp)

        gendrin_deg=wave.gendrin_angle(w[kk],wlhr,wce)


        ref_ind=np.sqrt(nx[kk]**2+ny[kk]**2+nz[kk]**2)

        Ytmp=np.absolute((np.rad2deg(psitmp)-theta_res_deg))


        psi[kk]=np.rad2deg(psitmp)
        theta_res[kk]=theta_res_deg
        Y[kk]=Ytmp
        S_stix[kk]=S_stixtmp
        D_stix[kk]=D_stixtmp
        P_stix[kk]=P_stixtmp
        R_stix[kk]=R_stixtmp
        L_stix[kk]=L_stixtmp
        gendrin_out[kk]=gendrin_deg
        wce_out[kk]=wce
        wcH_out[kk]=wcH
        wcHe_out[kk]=wcHe
        wcO_out[kk]=wcO
        wpe_out[kk]=wpe
        wpH_out[kk]=wpH
        wpHe_out[kk]=wpHe
        wpO_out[kk]=wpO    
        wlhr_out[kk]=wlhr
        wuhr_out[kk]=wuhr
        ref_out[kk]=ref_ind


    raytime = dt.datetime(2010,1,1,0,0,0)
    tmp_coords = coord.Coords(list(zip(posx, posy, posz)),'SM','car',units=['m','m','m'])
    tvec_datetime = [raytime + dt.timedelta(seconds=s) for s in time]
    tmp_coords.ticks = Ticktock(tvec_datetime)
    tmp_coords = tmp_coords.convert('MAG','sph')

    alt_ray=tmp_coords.radi
    lat_ray=tmp_coords.lati
    lon_ray=tmp_coords.long

    L_r = (alt_ray/env.const.Re)/pow(np.cos(np.deg2rad(lat_ray)),2)


    data={'time':time,
            'posx':posx,
            'posy':posy,
            'posz':posz,
            'vprelx':vprelx,
            'vprely':vprely,
            'vprelz':vprelz,
            'vgrelx':vgrelx,
            'vgrely':vgrely,
            'vgrelz':vgrelz,
            'nx':nx,
            'ny':ny,
            'nz':nz,
            'Bx':Bx,
            'By':By,
            'Bz':Bz,
            'w':w,
            'Ne':Ns1,
            'NH':Ns2,
            'NHe':Ns3,
            'NO':Ns4,
            'psi':psi,
            'theta_res':theta_res,
            'gendrin' : gendrin_out,
            'Y':Y,
            'L':L_r,
            'alt':alt_ray,
            'lat':lat_ray,
            'lon':lon_ray,
            'damp': damp_cut,
            'S_stix':S_stix,
            'D_stix':D_stix,
            'P_stix':P_stix,
            'R_stix':R_stix,
            'L_stix':L_stix,
            'wce':wce_out,
            'wcH':wcH_out,
            'wcHe':wcHe_out,
            'wcO':wcO_out,
            'wpe':wpe_out,
            'wpH':wpH_out,
            'wpHe':wpHe_out,
            'wpO':wpO_out, 
            'wlhr':wlhr_out,
            'wuhr':wuhr_out          
                }


    df = pd.DataFrame(data)

    df.to_csv(ray_file_name+'_appended.csv',sep=',',index=False)

    print('Saved file as:', ray_file_name+'_appended.csv')
