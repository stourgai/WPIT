import numpy as np

#####environment_mod.omega_cyclotron############################################

#Description:Routine to calculate the gyrofrequency of a particle
#Inputs:
# B_arg: Geomagnetic field strength in T
# q_arg: particle charge in Cb
# m_arg: particle mass in kg
#Outputs:
# omega_tmp: particle gyrofrequency in rad/s

#############################################################################


def omega_cyclotron(B_arg,q_arg,m_arg):
    #----calculate the cyclotron frequency
    #q_arg the species charge in Coulomb
    #m_arg the species mass in kg
    #B_arg the geomagnetic field strength in Tesla
    omega_tmp=(np.abs(q_arg)*B_arg)/m_arg
    return omega_tmp