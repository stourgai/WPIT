import numpy as np
import os
import sys
from scipy import integrate as scint
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)


import environment_mod as env
import waveproperties_mod as wave
import ray_mod as ray
import Landau_mod as landau

def hot_dispersion_imag(f, kperp, kpar, w, m, wch, qh, mh, qs, Ns, ms, nus, B0 ):

    eps=2.220446049250313e-16

    [S,D,P,R,L] = wave.stix_parameters(w, Ns[0],Ns[1],Ns[2],Ns[3], B0)
    

    integrand_vperp=lambda vperp: landau.integrand(f, vperp, kperp, kpar, w, m,
                                      wch, qh, mh,R, L, P, S)
    

    k=env.const.c_light
    integrand_vperpnorm= lambda vperp: k*integrand_vperp(k*vperp)

    integrand_t=lambda t: ((1+eps)/(t**2+eps))*integrand_vperpnorm((1-t+eps)/(t+eps))

    # integrated_integrand=scint.quad(integrand_t,0,1,epsabs=1e-3,maxp1=100,limit=100)
    integrated_integrand=scint.quad(integrand_vperpnorm,0,1,epsabs=1e-3,maxp1=100,limit=100)
    # integrated_integrand=scint.quad(integrand_vperp,0,k,maxp1=100,limit=100)
    ret=integrated_integrand[0]
  
    return ret

