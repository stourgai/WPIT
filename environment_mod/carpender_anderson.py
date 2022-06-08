import numpy as np
def L_pp(Kpmax):
    tmp=5.6-0.46*Kpmax
    return tmp

#saturated plasmasphere 2.25<=L<=Lppi
def sat_plasmasphere(L,d,R):
    log_ne=(-0.3145*L+3.9043)+(0.15*(np.cos(2*np.pi*(d+9)/365)-0.5*np.cos(4*np.pi*(d+9)/365))+0.00127*R-0.0635)*np.exp(-(L-2/1.5))
    ne=10**log_ne
#     print(ne)
    return ne

#Lppi<=L<=Lppo
def plasmapause(L,Lppi,ne_lppi,mlt,d,R):
    if mlt>=0 and mlt<6:
        ne=ne_lppi*10**(-(L-Lppi)/0.1)
    else:
        ne=ne_lppi*10**(-(L-Lppi)/(0.1+0.011*(mlt-6)))
    return ne

#extended plasma trough 2.25<=L<=8
def ext_trough(L,mlt):
    if mlt>=0 and mlt<6:
        ne=(5800+300*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    else:
        ne=(-800+1400*mlt)*L**(-4.5)+(1-np.exp(-(L-2)/10))
    return ne  

def trough(ne_Lppo,L,Lppo):
    ne=ne_Lppo*((L/Lppo)**(-4.5))+(1-np.exp(-(L-2)/10))
    return ne

#find the Lppo from the cross point of the two functions
def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc





def carpender_anderson(Lsh,Kpmax,day,mlt,Rb):

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
    # print(ne_Lppo)
#     for i in range(0,len(L_array)):
#         if L_array[i]<L_ppi:
#             ne=sat_plasmasphere(L_array[i],day,Rb)
#             ne_final.append(ne)
#             L_final.append(L_array[i])
#         if L_ppi<=L_array[i]<=L_ppo:
#             ne=plasmapause(L_array[i],L_ppi,ne_lppi,mlt,day,Rb)
#             ne_final.append(ne)
#             L_final.append(L_array[i])
#         if L_ppo<=L_array[i]<=8:
#             ne=trough(ne_Lppo,L_array[i],L_ppo)
#             ne_final.append(ne)
#             L_final.append(L_array[i]) 
    
    if Lsh<L_ppi:
        ne_eq=sat_plasmasphere(Lsh,day,Rb)
    if L_ppi<=Lsh<=L_ppo:
        ne_eq=plasmapause(Lsh,L_ppi,ne_lppi,mlt,day,Rb)
    if L_ppo<=Lsh<=8:
        ne_eq=trough(ne_Lppo,Lsh,L_ppo)
        
    return ne_eq