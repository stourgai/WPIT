def omega_plasma(n_arg,q_arg,m_arg):
    omegap_tmp=np.sqrt((n_arg*q_arg*q_arg)/(m_arg*const.epsilon0))
    fp_tmp=omegap_tmp/(2*np.pi)
    return omegap_tmp,fp_tmp
