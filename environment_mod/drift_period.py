import numpy as np
from environment_mod import const
def drift_period(B_arg,m_arg,v_arg,L_arg,aeq_arg):
    fac1=(2*np.pi*const.qe*B_arg*((const.Re)**3))/(m_arg*v_arg*v_arg)
    fac2=1/(L_arg*const.Re)
    fac3=1-(1/3)*((np.sin(aeq_arg))**0.62)
    tmp=fac1*fac2*fac3
    return tmp
