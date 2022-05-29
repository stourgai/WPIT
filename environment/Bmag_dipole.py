def Bmag_dipole(L_arg,lamda_arg):
    #----calculate the dipole magnetic field strength
    #L_arg geomagnetic L shell
    #lamda_arg geomagnetic latitude in rad
    
    B0=31200*10**(-9)
    slat=np.sin(lamda)
    clat=np.cos(lamda)
    slat_term = np.sqrt(1 + 3*slat*slat)
    Bmag=(B0/(L_shell**3))*slat_term/(clat**6)

    return Bmag
