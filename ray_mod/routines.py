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
        wps21 = (Ns[0]*qs[0]**2/ms[0]/env.const.eps0)
        wcs1 = ((qs[0]*B0mag)/ms[0])

        wps22 = (Ns[1]*qs[1]**2/ms[1]/env.const.eps0)
        wcs2 = ((qs[1]*B0mag)/ms[1])

        wps23 = (Ns[2]*qs[2]**2/ms[2]/env.const.eps0)
        wcs3 = ((qs[2]*B0mag)/ms[2])

        wps24 = (Ns[3]*qs[3]**2/ms[3]/env.const.eps0)
        wcs4 = ((qs[3]*B0mag)/ms[3])


        S_,D_,P_,R_,L_ = wave.stix_parameters(w, Ns[0],Ns[1],Ns[2],Ns[3], B0mag)

        tantheta_res2=-P_/S_
        tan_theta=np.sqrt(np.absolute(tantheta_res2))
        theta_res_tmp=np.arctan(tan_theta)

        theta_res_deg_tmp=np.rad2deg(theta_res_tmp)
        theta_res[ll]=theta_res_deg_tmp

        Ytmp=np.absolute((np.rad2deg(psi[ll])-theta_res_deg_tmp))

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

    L_r = (alt_ray/env.const.Re)/(np.cos(np.deg2rad(lat_ray))*np.cos(np.deg2rad(lat_ray)))

    ref_ind=np.sqrt(ray_in.nx**2+ray_in.ny**2+ray_in.nz**2)

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