import numpy as np
import os 
import sys
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from WPIT.Environment_mod import const
from  WPIT.LandauDamp_mod import distribution_bell,distribution_bortnik
#####Landau_mod.distribution_golden###############################################

#Description:Calculate the thermal electron distribution according to Golden et al.[2010]
#Inputs:
# vperp: perpendicular velocity
# vpar: parallel velocity
# scale: scale factor
#Outputs:
# f: electron distribution

#############################################################################

def distribution_golden(vperp,vpar,kpmax,Lmeas):

    fbell=distribution_bell(vperp,vpar)/(10**12)
    fbortnik=distribution_bortnik(vperp,vpar)/(10**12)

    Lpp=5.6-0.46*kpmax

    alpha=5
    w_bell=(np.exp(-alpha*(Lmeas-Lpp)))/(1+np.exp(-alpha*(Lmeas-Lpp)))
    w_bortnik=(np.exp(alpha*(Lmeas-Lpp)))/(1+np.exp(alpha*(Lmeas-Lpp)))


    e_fhybrid=(np.log(fbell)*w_bell+np.log(fbortnik)*w_bortnik)/(w_bortnik+w_bell)

    fhybrid=np.exp(e_fhybrid)*(10**12)

    return fhybrid
