import numpy as np


def res_angle(P_arg,S_arg):
    tansq=-P_arg/S_arg
    tan=np.sqrt(tansq)
    thetares=np.arctan(tan)
    return thetares

