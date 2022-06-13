import numpy as np
import os
import sys
from scipy import special as scp

current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import environment_mod as env
import waveproperties_mod as wave
import Landau_mod as landau

def integrand(f,vperp,kperp,kpar,w,m,wch,qh,mh,R,L,P,S):

    theta = np.arctan2(kperp,kpar)
    #refractive index
    n=np.sqrt((env.const.c_light**2/w**2)*(kperp**2+kpar**2))

    ct=np.cos(theta)
    st=np.sin(theta)

    vperp_vec=vperp
    ret=np.zeros(np.size(vperp_vec))

    


    # for vperp_i in range (0,len(vperp_vec)):

    for vperp_i in range (0,1):
        # vperp=vperp_vec[vperp_i]
        vperp=vperp_vec

        suminteg=0

        for mi in range(0,len(m)):
            Jm=scp.jv(m[mi], kperp*vperp/wch)
            Jm_m1=scp.jv(m[mi]-1, kperp*vperp/wch)
            Jm_p1=scp.jv(m[mi]+1, kperp*vperp/wch)
            

            Jm=np.real(Jm)
            Jm_m1=np.real(Jm_m1)
            Jm_p1=np.real(Jm_p1)

            #evaluate vpar only at resonance
            vpar = (w-m[mi]*wch)/kpar
            # evaluate G1
            G1 = landau.fG1(f, vperp, vpar, kperp, kpar, w)
            # evaluate G2
            G2 = landau.fG2(f, vperp, vpar, kperp, kpar, w, m[mi], wch)      

            suminteg=suminteg+  \
                (G1*((P-n**2*st**2)*(2*(L-n**2)*vperp*Jm_p1*Jm_p1+2*vperp*(R-n**2)*Jm_m1**2+    \
                    n**2*st*st*vperp*(Jm_p1-Jm_m1)*(Jm_p1-Jm_m1))-    \
                    n**2*ct*st*(2*vpar*Jm*(Jm_p1*(R-n**2)+Jm_m1*(L-n**2))+  \
                    n**2*ct*st*vperp*(Jm_p1-Jm_m1)*(Jm_p1-Jm_m1)))+  \
                    G2*(4*vpar*Jm*((L-n**2)*(R-n**2)+n*n*st*st*(S-n**2))-
                    2*n*n*ct*st*((R-n**2)*vperp*Jm_m1+(L-n**2)*vperp*Jm_p1)))
            
        ret[vperp_i]=-2*np.pi*np.pi*((qh*qh/mh/env.const.eps0)/(w*np.absolute(kpar)))*suminteg*vperp

    return ret