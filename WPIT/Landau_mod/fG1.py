import numpy as np
import os
import sys

current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import waveproperties_mod as wave

def fG1(f, vperp, vpar, kperp, kpar, w):

    # f = distribution function, a function handle that takes as input 
    #     (vperp,vpar) and returns the value of the distribution at that
    #     point.
    # vperp:  perpendicular velocity
    # vpar:   parallel velocity
    # kperp:  perpendicular component of k vector
    # kpar:   parallel component of k vector
    # w:      frequency, rad/s 

    eps=2.220446049250313e-16
    DEL=1e-8
    d=DEL*np.absolute(vperp)
    if d<10*eps:
        d=10*eps 

    df_dvperp = (f(vperp+d, vpar)-f(vperp-d,vpar))/(2*d)

    d=DEL*np.absolute(vpar)
    if d<10*eps:
        d=10*eps 

    df_dvpar = (f(vperp, vpar+d)-f(vperp,vpar-d))/(2*d)

    G1 = df_dvperp-(kpar/w)*(vpar*df_dvperp - vperp*df_dvpar)

    return G1