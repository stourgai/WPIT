"""
environment_mod.density_equ_carpenter_anderson

**Description**:
_____________________________________________________________________________________________________________________

Calculate the equatorial electron density according to Cartender and Anderson [1992] model

_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

Call model through environment_mod.density_equ_carpenter_anderson

**Reference**:

_____________________________________________________________________________________________________________________

Carpenter, D. L., and R. R. Anderson. "An ISEE/whistler model of equatorial electron density in the magnetosphere." 
Journal of Geophysical Research: Space Physics 97.A2 (1992): 1097-1108.

_____________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np





def L_pp(Kpmax):
    '''
    Auxiliary function for Carpenter-Anderson model

    Find plasmapause location
    '''
    tmp=5.6-0.46*Kpmax
    return tmp


def sat_plasmasphere(L,d,R):
    '''
    Auxiliary function for Carpenter-Anderson model

    Saturated plasmasphere 2.25<=L<=Lppi
    '''
    log_ne=(-0.3145*L+3.9043)+(0.15*(np.cos(2*np.pi*(d+9)/365)-0.5*np.cos(4*np.pi*(d+9)/365))+0.00127*R-0.0635)*np.exp(-(L-2/1.5))
    ne=10**log_ne

    return ne

def plasmapause(L,Lppi,ne_lppi,mlt,d,R):
    '''
    Auxiliary function for Carpenter-Anderson model

    Plasmapause segment Lppi<=L<=Lppo
    '''
    if mlt>=0 and mlt<6:
        ne=ne_lppi*10**(-(L-Lppi)/0.1)
    else:
        ne=ne_lppi*10**(-(L-Lppi)/(0.1+0.011*(mlt-6)))
    return ne


def ext_trough(L,mlt):
    '''
    Auxiliary function for Carpenter-Anderson model

    Extended plasma trough 2.25<=L<=8
    '''
    if mlt>=0 and mlt<6:
        ne=(5800+300*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    else:
        ne=(-800+1400*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    return ne  

def trough(ne_Lppo,L,Lppo):
    '''
    Auxiliary function for Carpenter-Anderson model

    Plasma trough Lppo<=L<=8
    '''
    ne=ne_Lppo*((L/Lppo)**(-4.5))+(1-np.exp(-(L-2)/10))
    return ne

#find the Lppo from the cross point of the two functions
def interpolated_intercept(x, y1, y2):
    """
    Routine to calculate the interceprion of two curves through interpolation
    """

    def intercept(p1, p2, p3, p4):
        """
        Find the intersection between two curves

        Needed fot the determination of the plasmapause outer limit Lppo
        as the intetcept point between the plasmapause segment and the
        extended plasma trough.


        **Outputs**: 
        
        The interception point (x,y) 
        """    

        def l(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersect_point(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = l([p1[0],p1[1]], [p2[0],p2[1]])
        L2 = l([p3[0],p3[1]], [p4[0],p4[1]])

        int_point = intersect_point(L1, L2)

        return int_point

    index = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xInt, yInt = intercept((x[index], y1[index]),((x[index+1], y1[index+1])), ((x[index], y2[index])), ((x[index+1], y2[index+1])))
    return xInt,yInt



def density_equ_carpenter_anderson(Lsh,Kpmax,day,mlt,Rb):
    """
    #####environment_mod.density_equ_carpenter_anderson##################################################################

    **Description**:
    _____________________________________________________________________________________________________________________

    Calculate the equatorial electron density according to Cartender and Anderson [1992] model

    _____________________________________________________________________________________________________________________
    _____________________________________________________________________________________________________________________

    **Inputs**:
    _____________________________________________________________________________________________________________________

    Lsh: L shell

    Kpmax: the maximum Kp value in the preceding 24 hours

    day: the day number

    mlt: the magnetic local time

    Rb: the 13-month average sunspot number

    _____________________________________________________________________________________________________________________
    _______________________________________________________________________________________________________________________

    **Outputs**:
    _____________________________________________________________________________________________________________________

    ne_eq: equatorial electron density in cm^-3

    _____________________________________________________________________________________________________________________
    ________________________________________________________________________________________________________________________

    **Reference**:

    _____________________________________________________________________________________________________________________

    Carpenter, D. L., and R. R. Anderson. "An ISEE/whistler model of equatorial electron density in the magnetosphere." 
    Journal of Geophysical Research: Space Physics 97.A2 (1992): 1097-1108.

    _____________________________________________________________________________________________________________________
    ________________________________________________________________________________________________________________________

    """


    L_plasma=[]
    ne_plasma=[]
    ne_trough=[]
    ne_final=[]
    L_final=[]
    # ne_final.append(10**6)
    # L_final.append(1)
    L_array=np.arange(2.25,8,0.001)

    L_ppi=L_pp(Kpmax)
    ne_lppi=sat_plasmasphere(L_ppi,day,Rb)

    for i in range(0,len(L_array)):
        if L_array[i]<L_ppi:
            ne=sat_plasmasphere(L_array[i],day,Rb)
            ne_plasma.append(ne)
            L_plasma.append(L_array[i])
        if L_array[i]>L_ppi:
            ne=plasmapause(L_array[i],L_ppi,ne_lppi,mlt,day,Rb)
            ne_plasma.append(ne)
            L_plasma.append(L_array[i])
    for i in range(0,len(L_array)):
        ne_tr=ext_trough(L_array[i],mlt)
        ne_trough.append(ne_tr)

    ne_plasma=np.asarray(ne_plasma)
    ne_trough=np.asarray(ne_trough)

    xc, yc = interpolated_intercept(L_array,ne_plasma,ne_trough)
    L_ppo=xc[0][0]

    ne_Lppo=plasmapause(L_ppo,L_ppi,ne_lppi,mlt,day,Rb)

    if Lsh<L_ppi:
        ne_eq=sat_plasmasphere(Lsh,day,Rb)
    if L_ppi<=Lsh<=L_ppo:
        ne_eq=plasmapause(Lsh,L_ppi,ne_lppi,mlt,day,Rb)
    if L_ppo<=Lsh<=8:
        ne_eq=trough(ne_Lppo,Lsh,L_ppo)
        
    return ne_eq