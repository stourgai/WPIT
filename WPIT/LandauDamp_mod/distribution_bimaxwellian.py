import numpy as np
import WPIT.Environment_mod as env
#####Landau_mod.bi_maxwellian_distribution###############################################

#Description:Calculate the anistorpic bi-maxwellian distribution
#Inputs:
# vperp: perpendicular velocity
# vpar: parallel velocity
#Outputs:
# f: electron distribution

#############################################################################

def distribution_bimaxwellian(vperp,vpar):

    nh=2*10**(-3) #hot electron density in m^-3


    Uthpar=0.05*env.const.c_light  #parallel component of thermal momentum per unit mass
    Uthper=0.03*env.const.c_light  #perpendicular component of thermal momentum per unit mass

    gamma=1/(np.sqrt(1-(vperp*vperp+vpar*vpar)/(env.const.c_light*env.const.c_light)))

    upar=gamma*vpar
    uper=gamma*vperp

    beta=0.01 #defines the loss cone, the larger the beta the larger the loss cone

    fac1=nh/(2*np.pi*(3/2)*Uthpar*Uthper*Uthper)
    fac2=np.exp(-(upar*upar)/(2*Uthpar*Uthpar))
    fac3=1/(1-beta)
    fac4=np.exp((-uper*uper)/(2*Uthper*Uthper))-np.exp((-uper*uper)/(2*beta*Uthper*Uthper))

    fmaxw=fac1*fac2*fac3*fac4


    return fmaxw
