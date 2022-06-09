import numpy as np

def densities_palsmasphere_sheeley(L_arg):
    ne_mean=1390*(3/L_arg)**4.8
    ne_min=ne_mean-440*(3/L_arg)**3.6
    ne_max=ne_mean+440*(3/L_arg)**3.6

    return ne_mean,ne_min,ne_max
