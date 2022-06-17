import numpy as np
import os
import sys

current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import waveproperties_mod as wave
import environment_mod as env
import Landau_mod as landau
def spatialdamping(f, kperp, kpar, w, m, wch, qh, mh,qs, Ns, ms, nus, B0):

# Compute the spatial damping rate ki in a hot plasma

    theta=np.arctan2(kperp,kpar)
    n=np.sqrt((env.const.c_light**2/w**2)*(kperp**2+kpar**2))

    [S,D,P,R,L] = wave.stix_parameters(w, Ns[0],Ns[1],Ns[2],Ns[3], B0)

    

    A = S*np.sin(theta)**2 + P*np.cos(theta)**2
    B = R*L*np.sin(theta)**2 + P*S*(1+np.cos(theta)**2)

    # print('A=',A)
    # print('B=',B)
    ki = 0
 
    for ii in range(0,len(f)):
        Di=landau.hot_dispersion_imag(f[ii], kperp, kpar, w, m , wch[ii],
                                    qh[ii], mh[ii], qs , Ns, ms, nus, B0)
        ki=ki + -(w/env.const.c_light)*(1/2)*(1/(4*n*(2*A*(n**2)-B)))*Di
        # ki=ki + -(w/const.c_light)*(1/(4*n*(2*A*(n**2)-B)))*Di #ST removed (1/2) from ki

    return ki