import numpy as np

import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
from WPIT.LandauDamp_mod.RayUtils_mod import read_appended_ray

def resonance_along_raypath(ray_file_name,mres,alpha_array):

    df=read_appended_ray(ray_file_name)

    time=df.time
    w=df.w
    nz=df.nz
    psi=df.psi
    wce=df.wce
    kpar=w*nz/env.const.c_light

    Ekev_res=np.zeros((len(alpha_array),len(time)),order='F')
    upar_res=np.zeros((len(alpha_array),len(time)),order='F')
    uper_res=np.zeros((len(alpha_array),len(time)),order='F')
    gamma_res=np.zeros((len(alpha_array),len(time)),order='F')
    Eres_res=np.zeros((len(alpha_array),len(time)),order='F')



    for k in range(0,len(alpha_array)):
        for i in range(0,len(time)):
            v_para_res, v_per_res, v_tot_res, E_res,gamm_res=wave.resonant_velocity(mres,w[i],kpar[i],wce[i],alpha_array[k],env.const.me)
            E_res = 5.105396765648739e5 *(1.0/np.sqrt( 1-(v_tot_res*v_tot_res/(env.const.c_light*env.const.c_light)) ) - 1 )/1000  
            upar_res[k,i]=v_para_res
            uper_res[k,i]=v_per_res
            gamma_res[k,i]=gamm_res
            Ekev_res[k,i]=E_res

    return time,psi,Ekev_res,upar_res,uper_res,gamma_res