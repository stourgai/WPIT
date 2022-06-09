import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const


def resonant_velocity(m_res_arg,w_wave_arg,kz_arg,wce_arg,alpha_arg,m_arg):
    fac1 = w_wave_arg*w_wave_arg*kz_arg*kz_arg
    fac2 = ((m_res_arg*wce_arg)*(m_res_arg*wce_arg))-(w_wave_arg*w_wave_arg)
    fac3 = kz_arg*kz_arg + ((m_res_arg*wce_arg)*(m_res_arg*wce_arg)) /((const.c_light*np.cos(alpha_arg))*(const.c_light*np.cos(alpha_arg)))

    if m_res_arg == 0:
        direction=-1.0*np.sign(kz_arg)
    else:
        direction = np.sign(kz_arg)*np.sign(m_res_arg)

    v_para_res = ( direction*np.sqrt(fac1 + fac2*fac3) - w_wave_arg*kz_arg) / fac3
    
    v_tot_res = v_para_res / np.cos(alpha_arg)
    v_per_res=v_tot_res*np.sin(alpha_arg)
    gamma_res=1.0/np.sqrt( 1-((v_tot_res*v_tot_res)/(const.c_light*const.c_light)))
    E_rest=m_arg*const.c_light*const.c_light

    E_res=(gamma_res-1)*E_rest
    return v_para_res, v_per_res, v_tot_res, E_res,gamma_res