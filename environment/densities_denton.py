def densities_denton(ne0_arg,lamda_arg):
    #-----calculate species densities (assuming 94%H+, 5.4%He+, 0.6%O+)
    #lamda geomagnetic latitude
    clat=np.cos(lamda_arg)
    n_e_tmp=ne0_arg*(clat**(-4))

    return n_e_tmp
