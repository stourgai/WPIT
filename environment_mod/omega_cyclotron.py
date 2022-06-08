import numpy as np
def omega_cyclotron(B_arg,q_arg,m_arg):
    #----calculate the cyclotron frequency
    #q_arg the species charge in Coulomb
    #m_arg the species mass in kg
    #B_arg the geomagnetic field strength in Tesla
    omega_tmp=(np.abs(q_arg)*B_arg)/m_arg
    return omega_tmp