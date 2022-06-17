import pandas as pd    
import numpy as np
import os
import sys
import datetime as dt
from spacepy import coordinates as coord
from spacepy.time import Ticktock

current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
sys.path.append(fpath)

import environment_mod as env
import waveproperties_mod as wave
import ray_mod as ray
import Landau_mod as landau

def read_appended_ray(ray_file_name):
    df=pd.read_csv(ray_file_name)

    

    return df