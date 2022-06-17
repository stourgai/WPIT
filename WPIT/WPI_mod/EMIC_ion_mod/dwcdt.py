import numpy as np

def dwcdt(ppar,m,gamma,dwcdz):
    tmp=(ppar/(m*gamma))*dwcdz
    return tmp