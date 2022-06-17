import pandas as pd    
import numpy as np

def read_input_ray(ray_file_name):


    ray_=pd.read_csv(ray_file_name,delim_whitespace=True, header=None)
    ray_.columns=['ray','ray_stop','time','posx','posy','posz','vprelx','vprely','vprelz',
                    'vgrelx','vgrely','vgrelz','nx','ny','nz','Bx','By','Bz','w','Nspec','qs1','qs2', 'qs3',
                    'qs4', 'ms1','ms2','ms3','ms4','Ns1','Ns2','Ns3','Ns4','nus1','nus2','nus3','nus4' ]

    timef=ray_.time
    posxf=ray_.posx
    posyf=ray_.posy
    poszf=ray_.posz
    vprelxf=ray_.vprelx
    vprelyf=ray_.vprely
    vprelzf=ray_.vprelz
    vgrelxf=ray_.vgrelx
    vgrelyf=ray_.vgrely
    vgrelzf=ray_.vgrelz
    nxf=ray_.nx
    nyf=ray_.ny
    nzf=ray_.nz
    Bxf=ray_.Bx
    Byf=ray_.By
    Bzf=ray_.Bz
    wf=ray_.w
    Nspecf=ray_.Nspec
    qs1f=ray_.qs1
    qs2f=ray_.qs2
    qs3f=ray_.qs3
    qs4f=ray_.qs4
    ms1f=ray_.ms1
    ms2f=ray_.ms2
    ms3f=ray_.ms3
    ms4f=ray_.ms4
    Ns1f=ray_.Ns1
    Ns2f=ray_.Ns2
    Ns3f=ray_.Ns3
    Ns4f=ray_.Ns4
    nus1f=ray_.nus1
    nus2f=ray_.nus2
    nus3f=ray_.nus3
    nus4f=ray_.nus4

    freq=wf/(2*np.pi)

    data={'time':timef,'posx':posxf,'posy':posyf,'posz':poszf,'vprelx':vprelxf,'vprely':vprelyf,'vprelz':vprelzf,
                   'vgrelx':vgrelxf,'vgrely':vgrelyf,'vgrelz':vgrelzf,'nx':nxf,'ny':nyf,'nz':nzf,'Bx':Bxf,'By':Byf,'Bz':Bzf,'w':wf,
                   'Nspec':Nspecf,'qs1':qs1f,'qs2':qs2f, 'qs3':qs3f,
                   'qs4':qs4f, 'ms1':ms1f,'ms2':ms2f,'ms3':ms3f,'ms4':ms4f,'Ns1':Ns1f,'Ns2':Ns2f,'Ns3':Ns3f,'Ns4':Ns4f,'nus1':nus1f,
                   'nus2':nus2f,'nus3':nus3f,'nus4':nus4f,'freq':freq}    
    df = pd.DataFrame(data)

    return df

