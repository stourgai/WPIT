import numpy as np

#####Landau_mod.bell_distribution###############################################

#Description:Calculate the electron distribution according to Bell [2002]
#Inputs:
# vperp: perpendicular velocity
# vpar: parallel velocity
#Outputs:
# f: electron distribution

#############################################################################

def distribution_bell(vperp,vpar):
    a=4.9e5
    b=8.3e14
    c=5.4e23

    v0=1

    v=100*np.sqrt(vperp*vperp+vpar*vpar+v0*v0)  #cm/s
    # v=100*(vm/(np.sqrt(1-(vm**2/env.const.c_light**2))))

    fbell=(a/(v**4))-(b/(v**5))+(c/(v**6))

    

    f=fbell*(10**12)  #convert to s^3/m^6 from s^3/cm^6

    f=10*f
    return f
