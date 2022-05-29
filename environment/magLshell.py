def magLshell(r_arg,lat_arg):
#     r=np.sqrt(posx[i]*posx[i]+posy[i]*posy[i]+posz[i]*posz[i])
    L_tmp=r_arg/(const.Re*np.cos(lat_arg)*np.cos(lat_arg))
    return L_tmp
