import numpy as np
from scipy.special import jn 

#####EMIC_ion_mod.dppardt###############################################

#Description:Routine to calculate the time derivative of the parallel momentum
#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# eta_arg: wave-particle phase in rad
# gamma_arg: Lorentz factor
# m_res_arg: Resonance order
# qi_arg: particle charge
# mi_arg: particle mass
# wtau_sq_arg: 𝜔2𝜏𝑚
# Bell parameter
# kz_arg: z component of the wave number
# Ewz_arg: z component of the electric field
# beta_arg: 𝛽
# wR_arg 𝜔𝑅
# wL_arg 𝜔𝐿
# Bmag_arg: magnitude of the geomagnetic field
# dBdz_arg: derivative of geomagnetic field with respect to the distance along the magnetic field line 

#Outputs:
# lrk: 𝑑𝑝∥𝑑𝑡
#############################################################################

def dppardt(pper_arg,eta_arg,gamma_arg,m_res_arg,qi_arg,mi_arg,Ewz_arg,beta_arg,wR_arg,wL_arg,Bmag_arg,dBdz_arg):
    tmpwave=((-1)**(m_res_arg+1))*(qi_arg*Ewz_arg*jn((m_res_arg),beta_arg)-wR_arg*pper_arg*jn((m_res_arg+1),beta_arg)+wL_arg*pper_arg*jn((m_res_arg-1),beta_arg))*np.sin(eta_arg)
    tmpadiabat=((pper_arg*pper_arg)/(2*gamma_arg*mi_arg*Bmag_arg))*dBdz_arg
    tmp=tmpwave-tmpadiabat
    return tmp