import numpy as np

#####EMIC_ion_mod.wpi_params###############################################

#Description:Routine to calculate needed parameters for EMIC_ion interactions
#Inputs:
# pper_arg: perpendicular to the geomagnetic field momentum
# kper_arg: perpendicular component of the wave number
# qi_arg: particle charge
# mi_arg: particle mass
# Bmag_arg: magnitude of the geomagnetic field
# Exw_arg: x component of the wave electric field in V/m
# Eyw_arg: y component of the wave electric field in V/m
# Bxw_arg: x component of the wave magnetic field in T
# Byw_arg: y component of the wave magnetic field in T
# gamma_arg: Lorentz factor

#Outputs:
# beta_tmp: 𝛽
# BwR: 𝐵𝑅𝑤
# BwL: 𝐵𝐿𝑤
# EwR: 𝐸𝑅𝑤
# EwL: 𝐸𝐿𝑤
# pwR: 𝑝𝑅𝑤
# pwL: 𝑝𝐿𝑤
# wR: 𝜔𝑅
# wL: 𝜔𝐿

#############################################################################

def wpi_params(pper_arg,kper_arg,qi_arg,mi_arg,Bmag_arg,Exw_arg,Eyw_arg,Bxw_arg,Byw_arg,gamma_arg):
    beta_tmp=-(kper_arg*pper_arg)/(qi_arg*Bmag_arg)
    BwR=0.5*(Bxw_arg+Byw_arg)
    BwL=0.5*(Bxw_arg-Byw_arg)
    EwR=0.5*(Exw_arg+Eyw_arg)
    EwL=0.5*(Exw_arg-Eyw_arg)    
    pwR=gamma_arg*mi_arg*(EwR/BwR)
    pwL=gamma_arg*mi_arg*(EwL/BwL)
    wR=(qi_arg*BwR)/(gamma_arg*mi_arg)
    wL=(qi_arg*BwL)/(gamma_arg*mi_arg)    
    
    return beta_tmp,BwR,BwL,EwR,EwL,pwR,pwL,wR,wL
