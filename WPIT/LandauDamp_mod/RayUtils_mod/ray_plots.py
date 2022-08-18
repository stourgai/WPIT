 
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from spacepy import coordinates as coord
from spacepy.time import Ticktock


import WPIT.Environment_mod as env
from WPIT.LandauDamp_mod.RayUtils_mod import read_appended_ray

def ray_plots(ray_file_name):

    df=read_appended_ray(ray_file_name)

    time=df.time
    posx=df.posx
    posy=df.posy
    posz=df.posz
    vphasex=df.vprelx
    vphasey=df.vprely
    vphasez=df.vprelz
    vgroupx=df.vgrelx
    vgroupy=df.vgrelx
    vgroupz=df.vgrelx
    nx=df.nx
    ny=df.ny
    nz=df.nz
    Bx=df.Bx
    By=df.By
    Bz=df.Bz
    w=df.w
    Ne=df.Ne
    NH=df.NH
    NHe=df.NHe
    NO=df.NO
    psi=df.psi
    theta_res=df.theta_res
    gendrin=df.gendrin
    Y=df.Y
    L=df.L
    alt=df.alt
    lat=df.lat
    lon=df.lon
    damp=df.damp
    S_stix=df.S_stix
    D_stix=df.D_stix
    P_stix=df.P_stix
    R_stix=df.R_stix
    L_stix=df.L_stix
    S_stix=df.S_stix
    wce=df.wce
    wcH=df.wcH
    wcHe=df.wcHe
    wcO=df.wcO
    wpe=df.wpe
    wpH=df.wpH
    wpHe=df.wpHe
    wpO=df.wpO
    wlhr=df.wlhr
    wuhr=df.wuhr

    ref_ind=np.sqrt(nx*nx+ny*ny+nz*nz)
    kpar=w*nz/env.const.c_light

    t_new=np.arange(0,max(time),0.0001)
    lat_int=np.interp(t_new,time,lat)
    L_int=np.interp(t_new,time,L)
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

    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('time (s)')
    ax1.set_ylabel('w ', color=color)
    ax1.set_title('Wave frequency, lower hybrid resonance frequency \n  and ray latitude as a function of time')
    ax1.plot(time,wlhr, color=color,label='$\omega_{LHR}$')
    ax1.plot(time,w, color="tab:pink",label='$\omega_{wave}$')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_xlim(0,np.max(time))
    ax1.grid(alpha=.3)
    ax1.legend()
    ax2 = ax1.twinx()  
    color = 'tab:blue'
    ax2.set_ylabel('lat(deg)', color=color) 
    ax2.plot(time,lat, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout() 
    ax1.legend()

    #----Power Plot------#
    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    ax2.set_ylim(0,1)
    ax2.set_xlim(0,np.max(time))
    ax2.set_title('Landau damping')
    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Normalised Wave Power')

    ax2.plot(time,damp,c='tab:blue', alpha=0.75)
    plt.grid(axis='both', alpha=.3)
    plt.show()

    #----L Plot------#
    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    # ax2.set_ylim(0,1)
    ax2.set_xlim(0,np.max(time))

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('L Shell')
    ax2.set_title('Ray L shell ')
    ax2.plot(time,L,c='tab:green', alpha=0.75)
    plt.grid(axis='both', alpha=.3)
    plt.show()

    #----Lat Plot------#
    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    # ax2.set_ylim(0,1)
    ax2.set_xlim(0,np.max(time))

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Latitude [deg]')
    ax2.set_title('Ray latitude ')
    ax2.plot(time,lat,c='tab:orange', alpha=0.75)
    plt.grid(axis='both', alpha=.3)
    plt.show()

    #----Psi Plot------#
    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    # ax2.set_ylim(0,1)
    ax2.set_xlim(0,np.max(time))

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('$\psi$ [deg]')
    ax2.set_title('Wave Normal Angle ')
    ax2.plot(time,psi,c='tab:green', alpha=0.75)
    plt.grid(axis='both', alpha=.3)
    plt.show()    

    #----Ref Plot------#
    fig2 = plt.figure()
    ax2=fig2.add_subplot(111)

    # ax2.set_ylim(0,1)
    ax2.set_xlim(0,np.max(time))

    ax2.set_xlabel('Time [sec]')
    ax2.set_ylabel('Refractive Index')
    ax2.set_title('Refractive Index')
    ax2.plot(time,ref_ind,c='tab:red', alpha=0.75)
    plt.grid(axis='both', alpha=.3)    

    #------Combined L-lat---#
    fig5 = plt.figure()
    ax=fig5.add_subplot(111)
    ax.set_xlabel("Time [s]",fontsize=14)
    ax.set_ylabel("L-shell",color='tab:green',fontsize=14)
    ax.set_title('Ray L shell and Latitude')
    ax.plot(time,L,color='tab:green')
    ax.grid(alpha=.3)
    ax.set_xlim(0,np.max(time))
    ax2=ax.twinx()
    # ax2.set_ylim(-20,20)
    ax2.plot(time,lat,color='tab:blue')
    ax2.set_ylabel("Latitude [deg]",color='tab:blue',fontsize=14)

    plt.show()

    #------Combined L-lat---#
    fig5 = plt.figure()
    ax=fig5.add_subplot(111)
    ax.set_xlabel("Time [s]",fontsize=14)
    ax.set_ylabel("Wave normal angle (deg)",color='tab:green',fontsize=14)
    ax.set_title('Ray Refractive Index and Latitude')
    ax.plot(time,psi,color='tab:green')
    ax.grid(alpha=.3)
    ax.set_xlim(0,np.max(time))
    ax2=ax.twinx()
    # ax2.set_ylim(-20,20)
    ax2.plot(time,lat,c='tab:red')
    ax2.set_ylabel("Latitude [deg]",c='tab:red',fontsize=14)

    plt.show()    

        #----------------Plot Ray path-----------

    
    D2R = (np.pi/180.0)
    R_E = 6371e3
    R_IONO=1000e3

    plotsize = 8   
    # L_shells = [1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,3, 4, 5, 6]    # Field lines to draw
    L_shells = [2,3, 4, 5, 6]
    fig, ax = plt.subplots(1,3)

    flashtime = dt.datetime(2010,1,1,0,0,0)
    tmp_coords = coord.Coords(list(zip(posx, posy, posz)),'SM','car',units=['m','m','m'])
    tvec_datetime = [flashtime + dt.timedelta(seconds=s) for s in time]
    tmp_coords.ticks = Ticktock(tvec_datetime)
    tmp_coords2=tmp_coords.convert('SM','sph')
    tmp_coords = tmp_coords.convert('MAG','car')

    ref_indexs=coord.Coords(list(zip(nx, ny, nz)),'SM','car')
    tvec_datetime2 = [flashtime + dt.timedelta(seconds=s) for s in time]
    ref_indexs.ticks = Ticktock(tvec_datetime2)
    ref_indexs2=ref_indexs.convert('SM','sph')
    ref_indexs=ref_indexs.convert('MAG','car')

    mag_field=coord.Coords(list(zip(Bx, By, Bz)),'SM','car')
    tvec_datetime3 = [flashtime + dt.timedelta(seconds=s) for s in time]
    mag_field.ticks = Ticktock(tvec_datetime3)
    mag_field=mag_field.convert('MAG','car')


    posx=tmp_coords.x
    posy=tmp_coords.y
    posz=tmp_coords.z
    nx=ref_indexs.x
    ny=ref_indexs.y
    nz=ref_indexs.z
    Bx=mag_field.x
    Bz=mag_field.y
    By=mag_field.z


    #plot a circle for the earth and the ionosphere
    for i in [0, 1, 2]:
        earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
        iono  = plt.Circle((0,0),(R_E + R_IONO)/R_E, color='c',alpha=0.5, zorder=99)
        ax[i].add_patch(earth)   
        ax[i].add_patch(iono)
        tmp_coords.sim_time = time

    # Plot the fieldlines (dipole model; could use something more complex)
    for Lf in L_shells:
        # Plot dipole field lines for both profile views
        lam = np.linspace(-80,80,181)
        L_r = Lf*pow(np.cos(lam*D2R),2)
        Lx  = L_r*np.cos(lam*D2R)
        Ly  = L_r*np.sin(lam*D2R)
        ax[1].plot(Lx,Ly,color='r',linewidth=1,linestyle='dashed')   # Field line
        ax[2].plot(Lx,Ly,color='r',linewidth=1,linestyle='dashed')   # Field line
        ax[2].plot(-Lx,Ly,color='r',linewidth=1,linestyle='dashed')  # Field line (other side)

        lam = np.linspace(-180,180,181)
        Lx2  = Lf*np.cos(lam*D2R)
        Ly2  = Lf*np.sin(lam*D2R)
        ax[0].plot(Lx2,Ly2,color='r',linewidth=1,linestyle='dashed')

    lw=2

    for i in range(0,360,10):
        tmp=i*np.pi/180
        qx=[np.cos(tmp),6*np.cos(tmp)]
        qy=[np.sin(tmp),6*np.sin(tmp)]
        ax[2].plot(qx,qy,'g--',linewidth=0.5)
        ax[1].plot(qx,qy,'g--',linewidth=0.5)

    nmagxz=np.sqrt(nx**2+nz**2)
    nmagyz=np.sqrt(ny**2+nz**2)
    nmag=np.sqrt(nx**2+ny**2+nz**2)
    Bmagxz=np.sqrt(Bx**2+Bz**2)
    Bmagyz=np.sqrt(By**2+Bz**2)
    Bmag=np.sqrt(Bx**2+By**2+Bz**2)

    ax[0].plot(posx/R_E, posy/R_E, linewidth=lw)
    if posy[0] < 0:
        ax[1].plot(posx/R_E, posz/R_E, linewidth=lw, zorder=101)
        for kk in range(0,len(posx),5):
            ax[1].quiver( posx[kk]/R_E,  posz[kk]/R_E, nx[kk]/nmagxz[kk], nz[kk]/nmagxz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)
            ax[2].quiver( posy[kk]/R_E,  posz[kk]/R_E, ny[kk]/nmagyz[kk], nz[kk]/nmagyz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)
    else:
        ax[1].plot(posx/R_E, posz/R_E, linewidth=lw, zorder=10)
        for kk in range(0,len(posx),5):
            ax[1].quiver(posx[kk]/R_E,  posz[kk]/R_E, nx[kk]/nmagxz[kk], nz[kk]/nmagxz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)
            # ax[1].quiver(posx[kk]/R_E,  posz[kk]/R_E, Bx[kk]/Bmagxz[kk], Bz[kk]/Bmagxz[kk], scale=25.,
            # width=0.003, headwidth=2., headlength=3., color='r')
            ax[2].quiver( posy[kk]/R_E,  posz[kk]/R_E, ny[kk]/nmagyz[kk], nz[kk]/nmagyz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)
            # ax[2].quiver( posy[kk]/R_E,  posz[kk]/R_E, By[kk]/Bmagyz[kk], Bz[kk]/Bmagyz[kk], scale=25.,
            # width=0.003, headwidth=2., headlength=3., color='r')
    if  posx[0] > 0:
        ax[2].plot(posy/R_E, posz/R_E, linewidth=lw, zorder=101)
    else:
        ax[2].plot(posy/R_E, posz/R_E, linewidth=lw, zorder=10)

    plt.title('Ray Path')
    ax[0].set_title('XY')
    ax[1].set_title('XZ')
    ax[2].set_title('YZ')
    ax[1].set_yticks([])
    ax[2].set_yticks([])

    ax[1].set_xlabel('L (R$_E$)')
    ax[0].set_ylabel('L (R$_E$)')

    ax[0].set_xlim([-plotsize, plotsize])
    ax[0].set_ylim([-plotsize, plotsize])
    ax[1].set_xlim([0, plotsize])
    ax[1].set_ylim([-plotsize/2, plotsize/2])
    ax[2].set_xlim([-plotsize, plotsize])
    ax[2].set_ylim([-plotsize, plotsize])

    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')
    ax[2].set_aspect('equal')
    
    plt.show()

    #-------------4x4 plot---------------------------------
    fig, axs = plt.subplots(2, 2)

    earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
    iono  = plt.Circle((0,0),(R_E + R_IONO)/R_E, color='c',alpha=0.5, zorder=99)
    axs[0,0].add_patch(earth)   
    axs[0,0].add_patch(iono)
    axs[0,0].set_xlim([-plotsize, plotsize])
    axs[0,0].set_ylim([-plotsize, plotsize])

    for Lf in L_shells:
        # Plot dipole field lines for both profile views
        lam = np.linspace(-80,80,181)
        L_r = Lf*pow(np.cos(lam*D2R),2)
        Lx  = L_r*np.cos(lam*D2R)
        Ly  = L_r*np.sin(lam*D2R)
        axs[0,0].plot(Lx,Ly,color='r',linewidth=1,linestyle='dashed')   # Field line
        axs[0,0].plot(-Lx,Ly,color='r',linewidth=1,linestyle='dashed')
    for i in range(0,360,10):
        tmp=i*np.pi/180
        qx=[np.cos(tmp),6*np.cos(tmp)]
        qy=[np.sin(tmp),6*np.sin(tmp)]
        #print(qx,qy)
        axs[0,0].plot(qx,qy,'g--',linewidth=0.5)

    if posy[0] < 0:
        sc=axs[0,0].scatter(posx/R_E, posz/R_E, s=4, c=damp, cmap='jet', zorder=101, edgecolor='none')


        for kk in range(0,len(posx),5):
            axs[0,0].quiver( posx[kk]/R_E,  posz[kk]/R_E, nx[kk]/nmagxz[kk], nz[kk]/nmagxz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)
    else:
        sc=axs[0,0].scatter(posx/R_E, posz/R_E,s=4, c=damp, cmap='jet',  zorder=10, edgecolor='none')
        for kk in range(0,len(posx),5):
            axs[0,0].quiver(posx[kk]/R_E,  posz[kk]/R_E, nx[kk]/nmagxz[kk], nz[kk]/nmagxz[kk], scale=25.,
            width=0.003, headwidth=2., headlength=3.)

    axs[0,0].scatter(posx[0]/R_E, posz[0]/R_E, s=60, c='yellow', zorder=101,edgecolor='black')
    plt.colorbar(sc,ax=axs[0,0])
    # plt.text(0.5, 0.5, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    axs[0,1].set_ylim(0,1)
    axs[0,1].set_xlim(0,np.max(time))

    axs[0,1].set_xlabel('Time [sec]')
    axs[0,1].set_ylabel('Normalised Wave Power')

    axs[1,1].set_xlabel('L-Shell')
    axs[1,1].set_ylabel('Enhancment Factor')

    axs[0,0].set_xlabel('L-shell')
    axs[0,0].set_ylabel('L-shell')

    axs[0,1].plot(time,damp,c='tab:green', alpha=0.75)
    axs[0,1].grid(axis='both', alpha=.3)


    axs[1,0].set_xlabel("Time [s]",fontsize=14)
    axs[1,0].set_ylabel("L-shell",color='tab:green',fontsize=14)
    axs[1,0].plot(time,L,color='tab:green')
    axs[1,0].grid(alpha=.3)
    axs[1,0].set_xlim(0,np.max(time))
    ax2=axs[1,0].twinx()
    ax2.plot(time,lat,color='tab:red')
    ax2.set_ylabel("Latitude [deg]",color='tab:red',fontsize=14)


    bins=np.arange(0.05,4.05,0.1)

    axs[1,1].hist(L_new,bins=bins,range=[0,4],density=True,weights=mag_new,color='tab:green')
    axs[1,1].grid()



    axs[0,0].set_xlim([0, plotsize-2])
    axs[0,0].set_ylim([-(plotsize-2)/2, (plotsize-2)/2])

    # axs[0,0].set_aspect('equal')
    # axs[0,1].set_aspect('equal')
    # axs[1,0].set_aspect('equal')
    # axs[1,1].set_aspect('equal')
    plt.tight_layout()
    plt.show()

    fig, axs = plt.subplots()

    earth = plt.Circle((0,0),1,color='0.5',alpha=1, zorder=100)
    iono  = plt.Circle((0,0),(R_E + R_IONO)/R_E, color='c',alpha=0.5, zorder=99)
    axs.add_patch(earth)   
    axs.add_patch(iono)
    axs.set_xlim([-plotsize, plotsize])
    axs.set_ylim([-plotsize, plotsize])

    for Lf in L_shells:
        # Plot dipole field lines for both profile views
        lam = np.linspace(-80,80,181)
        L_r = Lf*pow(np.cos(lam*D2R),2)
        Lx  = L_r*np.cos(lam*D2R)
        Ly  = L_r*np.sin(lam*D2R)
        axs.plot(Lx,Ly,color='r',linewidth=1,linestyle='dashed')   # Field line
        axs.plot(-Lx,Ly,color='r',linewidth=1,linestyle='dashed')
    for i in range(0,360,10):
        tmp=i*np.pi/180
        qx=[np.cos(tmp),6*np.cos(tmp)]
        qy=[np.sin(tmp),6*np.sin(tmp)]
        #print(qx,qy)
        axs.plot(qx,qy,'g--',linewidth=0.5)

    if posy[0] < 0:
        sc=axs.scatter(posx/R_E, posz/R_E, s=4, c=damp, cmap='jet', zorder=101, edgecolor='none')


        for kk in range(0,len(posx),20):
            axs.quiver( posx[0]/R_E,  posz[0]/R_E, nx[0]/nmagxz[0], nz[0]/nmagxz[0], scale=15.,
            width=0.003, headwidth=2., headlength=3.)
    else:
        sc=axs.scatter(posx/R_E, posz/R_E,s=4, c=damp, cmap='jet',  zorder=10, edgecolor='none')
        for kk in range(0,len(posx),20):
            axs.quiver(posx[0]/R_E,  posz[0]/R_E, nx[0]/nmagxz[0], nz[0]/nmagxz[0], scale=15.,
            width=0.003, headwidth=2., headlength=3.)

    axs.scatter(posx[0]/R_E, posz[0]/R_E, s=60, c='yellow', zorder=101,edgecolor='black')
    plt.colorbar(sc,ax=axs)
    # plt.text(0.5, 0.5, 'matplotlib', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)


    axs.set_xlabel('L-shell')
    axs.set_ylabel('L-shell')
    axs.set_title('Ray Path \n with color coded Landau Damping')
#     axs.set_xlim([0, plotsize])
#     axs.set_ylim([-(plotsize)/2, (plotsize)/2])
    axs.set_xlim([2.5, 5.1])
    axs.set_ylim([-2, 2])
    # axs[0,0].set_aspect('equal')
    # axs[0,1].set_aspect('equal')
    # axs[1,0].set_aspect('equal')
    # axs[1,1].set_aspect('equal')
    plt.tight_layout()
    plt.show()