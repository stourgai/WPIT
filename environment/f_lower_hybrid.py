def f_lower_hybrid(fc_arg,fp_arg,ne_arg,nH_arg,nO_arg,nHe_arg):
    ion_fac=(((nH_arg)/const.mH)+((nHe_arg)/const.mHe)+((nO_arg)/const.mO))
    m_fac=(const.me/(ne_arg))
    tmp=np.sqrt(m_fac*((fc_arg*fc_arg*fp_arg*fp_arg)/(fc_arg*fc_arg+fp_arg*fp_arg)))
    return tmp
