import numpy as np
import os 
import sys
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

from WPIT.Environment_mod import const
#####Landau_mod.bortnik_distribution###############################################

#Description:Calculate the thermal electron distribution according to Bortnik et al.[2007,a]
#Inputs:
# vperp: perpendicular velocity
# vpar: parallel velocity
# scale: scale factor
#Outputs:
# f: electron distribution

#############################################################################

def distribution_bortnik(vperp,vpar):
    ko=6.25e11
    mdot=ko*const.me

    a1=0.755
    a0=np.log10(2.14*(10**7))

    nu=2*a1+2

    An=(2*10**(a0))/((0.5*mdot)**(a1-1))


    v=100*np.sqrt(vperp*vperp+vpar*vpar)  #cm/s


    fbort=An/(v**nu)



    f=fbort*(10**12)  #convert to s^3/m^6 from s^3/cm^6

    return f
