import os
import sys
import numpy as np
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from environment_mod import const

#####waveproperties_mod.ref_index##############################

#Description:Routine to calculate the the refractive index and the wave numbers
#Inputs:
# theta: wave normal angle in rad
# S:Stix S parameter
# P:Stix P parameter
# R:Stix R parameter
# L:Stix L parameter
# D:Stix D parameter
# w_wave_arg: wave frequency
#Outputs:
# eta_sq_plus: the plus (+) root of the dispersion relation
# eta_sq_minus: the minus (-) root of the dispersion relation
# ref_ind: refractive index (for the 𝑛2>0 root)
# kappa: wave number
# kappa_par: parallel component of the wave number
# kappa_per: perpendicular component of the wave number

##############################################################

def ref_index(theta_arg,w_wave_arg,S_arg,P_arg,R_arg,L_arg):
    A=S_arg*np.sin(theta_arg)*np.sin(theta_arg)+P_arg*np.cos(theta_arg)*np.cos(theta_arg)
    B=R_arg*L_arg*np.sin(theta_arg)*np.sin(theta_arg)+P_arg*S_arg*(1+np.cos(theta_arg)*np.cos(theta_arg))
    C=P_arg*R_arg*L_arg
    F=np.sqrt(B*B-4*A*C)
    eta_sq_plus=(B+F)/(2*A)
    eta_sq_minus=(B-F)/(2*A)
    if eta_sq_plus>0:
        ref_ind=np.sqrt(eta_sq_plus)
    else:
        ref_ind=np.sqrt(eta_sq_minus)

    
    kappa=ref_ind*w_wave_arg/const.c_light
    kappa_par=kappa*np.cos(theta_arg)
    kappa_per=kappa*np.sin(theta_arg)
    return eta_sq_plus,eta_sq_minus,ref_ind,kappa,kappa_par,kappa_per