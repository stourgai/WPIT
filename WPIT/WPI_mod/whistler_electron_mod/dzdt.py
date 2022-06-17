import numpy as np

def dzdt(ppar_arg,gamma_arg,mi_arg):
    krk=ppar_arg/(gamma_arg*mi_arg)
    return krk
