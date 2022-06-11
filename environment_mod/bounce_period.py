import numpy as np
from environment_mod import const

#####environment_mod.bounce_period###############################################

#Description:**Calculate the bounce period of a trapped particle**
#Inputs:
# **L_arg**: L shell
# **v_arg**: particle velocity in m/s
# **aeq_arg**: equatorial pitch angle
#Outputs:
# **tmp**: Bounce period in s

#################################################################################

def bounce_period(L_arg,v_arg,aeq_arg):
    fac1=0.117*L_arg*const.c_light/v_arg
    fac2=1-0.4635*((np.sin(aeq_arg))**(3/4))
    tmp=fac1*fac2
    return tmp
