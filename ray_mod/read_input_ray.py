import pandas as pd    
import numpy as np

def read_input_ray(ray_file_name):


    ray_in=pd.read_csv(ray_file_name,delim_whitespace=True, header=None)
    ray_in.columns=['ray','ray_stop','time','posx','posy','posz','vprelx','vprely','vprelz',
                    'vgrelx','vgrely','vgrelz','nx','ny','nz','Bx','By','Bz','w','Nspec','qs1','qs2', 'qs3',
                    'qs4', 'ms1','ms2','ms3','ms4','Ns1','Ns2','Ns3','Ns4','nus1','nus2','nus3','nus4' ]

    timef=ray_in.time
    posxf=ray_in.posx
    posyf=ray_in.posy
    poszf=ray_in.posz
    vprelxf=ray_in.vprelx
    vprelyf=ray_in.vprely
    vprelzf=ray_in.vprelz
    vgrelxf=ray_in.vgrelx
    vgrelyf=ray_in.vgrely
    vgrelzf=ray_in.vgrelz
    nxf=ray_in.nx
    nyf=ray_in.ny
    nzf=ray_in.nz
    Bxf=ray_in.Bx
    Byf=ray_in.By
    Bzf=ray_in.Bz
    wf=ray_in.w
    Nspecf=ray_in.Nspec
    qs1f=ray_in.qs1
    qs2f=ray_in.qs2
    qs3f=ray_in.qs3
    qs4f=ray_in.qs4
    ms1f=ray_in.ms1
    ms2f=ray_in.ms2
    ms3f=ray_in.ms3
    ms4f=ray_in.ms4
    Ns1f=ray_in.Ns1
    Ns2f=ray_in.Ns2
    Ns3f=ray_in.Ns3
    Ns4f=ray_in.Ns4
    nus1f=ray_in.nus1
    nus2f=ray_in.nus2
    nus3f=ray_in.nus3
    nus4f=ray_in.nus4

    freq=wf/(2*np.pi)

    data={'time':timef,'posx':posxf,'posy':posyf,'posz':poszf,'vprelx':vprelxf,'vprely':vprelyf,'vprelz':vprelzf,
                   'vgrelx':vgrelxf,'vgrely':vgrelyf,'vgrelz':vgrelzf,'nx':nxf,'ny':nyf,'nz':nzf,'Bx':Bxf,'By':Byf,'Bz':Bzf,'w':wf,
                   'Nspec':Nspecf,'qs1':qs1f,'qs2':qs2f, 'qs3':qs3f,
                   'qs4':qs4f, 'ms1':ms1f,'ms2':ms2f,'ms3':ms3f,'ms4':ms4f,'Ns1':Ns1f,'Ns2':Ns2f,'Ns3':Ns3f,'Ns4':Ns4f,'nus1':nus1f,
                   'nus2':nus2f,'nus3':nus3f,'nus4':nus4f,'freq':freq}    
    df = pd.DataFrame(data)

    return df

