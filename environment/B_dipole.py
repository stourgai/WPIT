def B_dipole(L_arg,lamda_arg,r_arg):
    import const
    B0=31200*10**(-9)
    Br_tmp=-2*B0*((const.Re/r_arg)**3)*np.sin(lamda_arg)
    Bl_tmp=B0*((const.Re/r_arg)**3)*np.cos(lamda_arg)
    Bt_tmp=0
    Bmag_tmp=np.sqrt(Br_tmp*Br_tmp+Bl_tmp*Bl_tmp+Bt_tmp*Bt_tmp)
    return Br_tmp,Bl_tmp,Bt_tmp,Bmag_tmp
