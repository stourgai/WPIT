import numpy as np
import os
import sys
from scipy import special as scp
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import waveproperties_mod as wave

def fG2(f, vperp, vpar, kperp, kpar, w, m, wch):

    # f = distribution function, a function handle that takes as input 
    #     (vperp,vpar) and returns the value of the distribution at that
    #     point.
    # vperp:  perpendicular velocity
    # vpar:   parallel velocity
    # kperp:  perpendicular component of k vector
    # kpar:   parallel component of k vector
    # w:      frequency, rad/s 
    # m:      resonance (0=landau, -1,+1=cyclotron, etc.)
    # wch:    hot (signed) gyrofrequency, negative for electrons
    DEL=1e-8

    eps=2.220446049250313e-16
    d=DEL*np.absolute(vperp)
    if d<10*eps:
        d=10*eps 

    df_dvperp = (f(vperp+d, vpar)-f(vperp-d,vpar))/(2*d)

    d=DEL*np.absolute(vpar)
    if d<10*eps:
        d=10*eps 

    df_dvpar = (f(vperp, vpar+d)-f(vperp,vpar-d))/(2*d)

    Jm=np.real(scp.jv(m,kperp*vperp/wch))

    G2 = Jm*(df_dvpar-(m*wch+eps)/(w*vperp+eps)*(vpar*df_dvperp - vperp*df_dvpar))
    return G2
