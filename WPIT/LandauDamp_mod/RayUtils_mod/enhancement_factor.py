
import numpy as np

import matplotlib.pyplot as plt

import WPIT.Environment_mod as env
from .read_appended_ray import read_appended_ray

def enhancement_factor(ray_file_name):

    df=read_appended_ray(ray_file_name)

    time=df.time
    w=df.w
    nz=df.nz
    psi=df.psi
    wce=df.wce
    kpar=w*nz/env.const.c_light
    lat=df.lat
    damp=df.damp
    Lf=df.L

    t_new=np.arange(0,max(time),0.0001)
    lat_int=np.interp(t_new,time,lat)
    L_int=np.interp(t_new,time,Lf)
    mag_int=np.interp(t_new,time,damp)

    time_new=[]
    lat_new=[]
    L_new=[]
    mag_new=[]
    time_new.append(t_new[0])
    lat_new.append(lat_int[0])
    L_new.append(L_int[0])
    mag_new.append(mag_int[0])
    for i in range(1,len(t_new)-1):
        if lat_int[i]>0.0 and lat_int[i+1]<0.0 or lat_int[i]<0.0 and lat_int[i+1]>0.0 :
            time_new.append(t_new[i])
            lat_new.append(lat_int[i])
            L_new.append(L_int[i])
            mag_new.append(mag_int[i])


        #-----LBin Plot-----#
    fig=plt.figure()
    fig.add_subplot(111)


    bins=np.arange(np.mean(L_new)-2,np.mean(L_new)+2,0.1)

    plt.hist(L_new,bins=bins,range=[0,4],density=True,weights=mag_new,color='tab:blue')
    # sns.displot(L_new,weights=mag_new,kde=False, discrete=True, bins=bins)
    plt.xlabel('L shell')
    plt.ylabel('Cavity enhancment')
    plt.grid()

    plt.show()
