import numpy as np

#####parallel_EMIC.nonlinear_S###############################################

#Inputs:
# H: non linear parameter H
# wtsq: trapping frequency squared

#Outputs:
# tmp: s

#############################################################################

def nonlinear_S(H_arg,wtsq_arg):
    tmp=H_arg/wtsq_arg
    return tmp